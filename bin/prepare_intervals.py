#!/usr/bin/env python3

# Usage: python prepare_intervals.py --bed /home/pmoris/itg/datasets/reference-genomes/Pfalciparum/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.bed --min-contig-size 2000000 --target-bin-size 2000000 --orphan-fraction .3

"""
Bin genomic regions/intervals/contigs into approximately balanced groups for parallel processing.

Important options:

- Standalone contigs:
    Contigs with size >= --min-contig-size/--standalone-contig-size are considered too large
    to be grouped with others and are therefore placed into their own bin.
    Defaults to total genome size / 20 (to keep the number of parallel bins reasonable).

- Groupable contigs:
    Contigs smaller than --min-contig-size are grouped using a First-Fit Decreasing (FFD)
    bin-packing heuristic, subject to a hard upper limit (--target-bin-size) on the total size
    of the bin. Defaults to size of the largest contig.

- --orphan-fraction:
    Very small bins (fraction of target-bin-size) may be merged to avoid left-over tiny bin.

--max-bins:
    Optional upper cap on the total number of bins, enforcing additional merging.

"""

import sys
import argparse
import logging
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Bin genomic intervals using first-fit decreasing for parallel processing"
    )

    parser.add_argument(
        "--bed",
        required=True,
        type=Path,
        help="Input BED file with genomic intervals",
    )
    parser.add_argument(
        "--min-contig-size",
        "--standalone-contig-size",
        dest="min_contig_size",
        type=int,
        default=None,
        help="Contigs >= this size (bp) get their own bin. Defaults to genome size / 20.",
    )
    parser.add_argument(
        "--target-bin-size",
        type=int,
        # required=True,
        default=None,
        help="Target bin size (bp) for grouping small contigs (= maximum size of bins that may contain multiple contigs). Defaults to size of the largest contig.",
    )
    parser.add_argument(
        "--orphan-fraction",
        type=float,
        default=0.2,
        help="Merge bins smaller than this fraction of target size (default: 0.2)",
    )
    parser.add_argument(
        "--max-bins",
        type=int,
        default=None,
        help="Maximum number of bins to produce (default: no limit)",
    )
    parser.add_argument(
        "--sort-mode",
        type=str,
        default="size",
        help="How to sort contigs inside each bin: 'size' or 'genomic' (default: 'size')",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO)",
    )

    return parser.parse_args()


def read_contigs(bed_file: Path, logger: logging.Logger):
    contigs = []

    with bed_file.open() as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            chrom, start, end = line.rstrip().split("\t")[:3]
            contigs.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "size": int(end) - int(start),
                }
            )

    logger.debug(f"Read {len(contigs)} contigs from {bed_file}")

    return contigs


def new_bin(contigs=None):
    contigs = contigs or []
    return {
        "contigs": contigs,
        "size": sum(c["size"] for c in contigs),
    }


def first_fit_decreasing(contigs, bins, target_size, logger):
    for contig in contigs:
        placed = False
        for b in bins:
            # Do not add small contigs to bins already exceeding target size
            if b["size"] > target_size:
                continue

            if b["size"] + contig["size"] <= target_size:
                b["contigs"].append(contig)
                b["size"] += contig["size"]
                placed = True
                break
        else:
            bins.append(new_bin([contig]))
            logger.debug(
                f"Created new bin for contig {contig['chrom']} ({contig['size']:,} bp)"
            )


def merge_orphan_bin(bins, target_size, orphan_fraction, logger):
    if len(bins) < 2:
        return

    bins.sort(key=lambda b: b["size"])
    smallest = bins[0]

    threshold = target_size * orphan_fraction
    if smallest["size"] < threshold:
        logger.info(
            f"Merging orphan bin ({smallest['size']:,} bp) into next smallest bin"
        )
        bins[1]["contigs"].extend(smallest["contigs"])
        bins[1]["size"] += smallest["size"]
        bins.pop(0)
    else:
        logger.warning(
            f"Smallest bin ({smallest['size']:,} bp) too large to be merged with next smallest bin"
        )


def enforce_max_bins(bins, max_bins, logger):
    """
    Reduce the number of bins to max_bins by repeatedly merging
    the two smallest bins.
    """
    if max_bins is None or len(bins) <= max_bins:
        return

    logger.info(
        f"Reducing bin count from {len(bins)} to {max_bins} by merging smallest bins"
    )

    # Work on a copy-like structure but mutate in place
    while len(bins) > max_bins:
        bins.sort(key=lambda b: b["size"])
        b1 = bins.pop(0)
        b2 = bins.pop(0)
        merged = {
            "contigs": b1["contigs"] + b2["contigs"],
            "size": b1["size"] + b2["size"],
        }
        bins.append(merged)
        logger.debug(
            f"Merged bins of size {b1['size']:,} bp and {b2['size']:,} bp "
            f"→ {merged['size']:,} bp"
        )


def write_outputs(bins, logger):
    logger.info(f"Creating {len(bins)} interval files")

    stats_lines = ["Bin\tNum_Contigs\tTotal_Size\tContigs\n"]

    for i, b in enumerate(bins, start=1):
        filename = f"interval_{i:03d}.bed"
        with open(filename, "w") as out:
            for c in b["contigs"]:
                out.write(f"{c['chrom']}\t{c['start']}\t{c['end']}\n")

        contig_names = ", ".join(c["chrom"] for c in b["contigs"])
        logger.info(f"Bin {i:3d}: {len(b['contigs']):2d} contigs, {b['size']:10,} bp")

        stats_lines.append(f"{i}\t{len(b['contigs'])}\t{b['size']}\t{contig_names}\n")

    with open("binning_stats.txt", "w") as stats:
        stats.writelines(stats_lines)


def validate_args(args, contigs, logger):
    sizes = [c["size"] for c in contigs]
    min_contig = min(sizes)

    if args.target_bin_size < min_contig:
        logger.warning(
            f"target-bin-size ({args.target_bin_size}) is smaller than the smallest contig ({min_contig})!"
            "All contigs will be forced into standalone bins...",
        )

    if args.min_contig_size > args.target_bin_size:
        logger.warning(
            f"min-contig-size ({args.min_contig_size}) is larger than target-bin-size ({args.target_bin_size})! "
            "All grouped bins will be smaller than any standalone contig..."
        )

    if not (0.0 < args.orphan_fraction < 1.0):
        logger.error("--orphan-fraction must be between 0 and 1")
        sys.exit(1)

    if args.max_bins is not None and args.max_bins <= 0:
        logger.error("--max-bins must be a positive integer")
        sys.exit(1)


def main():
    args = parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
    )
    logger = logging.getLogger(__name__)

    contigs = read_contigs(args.bed, logger)
    if not contigs:
        logger.error("No valid intervals found in BED file")
        sys.exit(1)

    # set default minimum contig size if not provided
    genome_size = sum(c["size"] for c in contigs)
    if args.min_contig_size is None:
        args.min_contig_size = genome_size // 20
        logger.info(
            f"min-contig-size not specified; using {args.min_contig_size} bp (~genome_size / 20)"
        )

    # set target bin size to largest contig if not provided
    if args.target_bin_size is None:
        args.target_bin_size = max(c["size"] for c in contigs)
        logger.info(
            f"target-bin-size was not specified; using largest contig size ({args.target_bin_size} bp)"
        )

    validate_args(args, contigs, logger)

    large_contigs = [c for c in contigs if c["size"] >= args.min_contig_size]
    small_contigs = [c for c in contigs if c["size"] < args.min_contig_size]

    logger.info(f"Total contigs: {len(contigs)}")
    logger.info(f"Minimum requested contig size: {args.min_contig_size:,} bp")
    logger.info(f"Target bin size: {args.target_bin_size:,} bp")
    logger.info(f"Large contigs: {len(large_contigs)}")
    logger.info(f"Small contigs: {len(small_contigs)}")

    # Initialize bins with large contigs - these do not need to be merged
    bins = [new_bin([c]) for c in large_contigs]

    # Sort small contigs by decreasing size for first-fit-decreasing bin packing
    small_contigs.sort(key=lambda c: c["size"], reverse=True)

    first_fit_decreasing(small_contigs, bins, args.target_bin_size, logger)

    # Merge orphan bins
    merge_orphan_bin(
        bins,
        args.target_bin_size,
        args.orphan_fraction,
        logger,
    )

    # Enforce maximum bin count
    enforce_max_bins(bins, args.max_bins, logger)

    # sanity checks
    all_assigned = [
        (c["chrom"], c["start"], c["end"]) for b in bins for c in b["contigs"]
    ]

    if len(all_assigned) != len(set(all_assigned)):
        logger.error("Duplicate contigs detected across bins!")
        sys.exit(1)

    if len(all_assigned) != len(contigs):
        logger.error(f"Lost contigs: expected {len(contigs)}, got {len(all_assigned)}")
        sys.exit(1)

    # Order contigs within bins
    for b in bins:
        b["contigs"].sort(key=lambda c: (c["chrom"], c["start"], c["end"]))

    # Order bins according to "genomic position" (=alphabetical) or size
    if args.sort_mode == "genomic":
        bins.sort(
            key=lambda b: (
                b["contigs"][0]["chrom"],
                b["contigs"][0]["start"],
                b["contigs"][0]["end"],
            )
        )
    elif args.sort_mode == "size":
        bins.sort(key=lambda b: b["size"], reverse=True)
    else:
        logger.error("Incorrect sort mode was selected.")
        sys.exit(1)

    sizes = [b["size"] for b in bins]
    logger.info("-" * 60)
    logger.info(f"Total bins: {len(bins)}")
    logger.info(f"Min bin size: {min(sizes):,} bp")
    logger.info(f"Max bin size: {max(sizes):,} bp")
    logger.info(f"Mean bin size: {sum(sizes) // len(sizes):,} bp")
    logger.info(f"Total genome: {sum(sizes):,} bp")

    write_outputs(bins, logger)


if __name__ == "__main__":
    main()
