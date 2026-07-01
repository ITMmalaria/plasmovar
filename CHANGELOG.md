# ITMmalaria/plasmovar: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html). 
Additionally, it (mostly) follows the [nf-core recommendations](https://nf-co.re/docs/specifications/pipelines/requirements/semantic_versioning).

## 0.4.0 - [date]

### `Added`

- Add option to emit invariant sites during joint genotype calling (`gatk GenotypeGVCFs -all-sites`) and enable it by default.

### `Improved`

- Changed default interval_list scatter count from 20 to 10 (which is what e.g. _P. falciparum_ 3D7 defaults to for its 14+2 chromosomes when using the `BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW` option for `IntervalListTools`).
- Run VQSR on sites-only VCF to speed up model building step (see [GATK VQSR tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering)).
- Improved documentation and citation message (e.g., clarify that this pipeline utilises nf-core components, but that it is not an official nf-core pipeline itself).
- Attempt to make Nextflow caching behaviour more deterministic by sorting the content of various channels before aggregating/combining them. E.g., scattered per-interval files are sorted on interval_list index prior to gathering for certain `GATK` steps, and samples are sorted prior to joining by intervals as the input for `GenomicsDBImport`.
- Removed `GATK4_HAPLOTYPECALLER` ext.args that was the same as the `HaplotypeCaller`'s default value (`--contamination-fraction-to-filter 0.0`).

### `Fixed`

- Added missing `an MQRankSum` option for VQSR INDEL recalibration model building.
- Reorder VQSR steps, clarify channel names and fix inputs for `VariantRecalibrator`/`ApplyVQSR`. Now both the SNP and INDEL model building steps are executed prior to recalibration. The INDEL recalibration model building step was using the SNP-recalibrated VCF as an input, instead of the raw VCF. (This was not wrong, but it was confusing and differed from the procedure outlined in the [GATK VQSR tutorial](https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering).)

## 0.3.0 - [2026-06-15]

### `Added`

- Allow order of host read removal and trimming to be toggled via `--trimming_before_hostremoval`.
- Add bgzip and tabix step after SnpEff annotation.
- Base quality score recalibration (BQSR), or more specifically the `gatk BaseRecalibrator` step, can now be run without scattering/parallellization over intervals (`--bqsr_scatter false`), or by scattering over distinct intervals from those used in the rest of the pipeline (`--bqsr_bed <bed-file>`). This avoids problems when the known sites VCF contains different genomic regions from those used elsewhere in the pipeline.
- Added option to save BQSR BAM files (`--save_bqsr_bam`).

### `Improved`

- Renamed various parameters:
   - `--save_intermediate_gvcf` -> `--save_gvcf`
   - `--save_bam_final` -> `--save_bam`
- Changed file publishing behaviour for:
   - VCF files: only the final bgzipped-filter-added-and-annotated VCF file is now published by default. All other intermediate VCF files (GenotypeGVCFs, MergeVcfs, filtered VCFs) are now only published when the `--save_intermediate_vcf` option is enabled.
   - VQSR output folder is now named `vqsr` instad of `filter_vqsr`.
- Switch resource requirements to VariantRecalibrator to low.

## v0.2.0 - [2026-06-10]

### `Improved`

- Improve lane and flowcell detection and validation during samplesheet parsing.
   - New boolean parameters `--strict_lane_detection` and `--strict_flowcell_detection` were introduced to replace the previous `--auto_detect_*` parameters. Detection is now always performed, but in strict mode, the pipeline will halt rather than merely show warnings when there are mismatches between the samplesheet-provided and detected values, or if no values were provided and detection fails. Mismatches between the detected values of read file pairs always result in an error.
   - Lane detection in fastq file name now returns the exact lane info, without number padding, and including any `L`-type prefix. This avoids false positive mismatch reporting between detected and provided values.
   - Fallback lane detection in fastq header instead of file name was added. This returns the actual value detected in the header (usually a single digit), cast as a string, instead of `L###`, which hopefully results in fewer false positive mismatch warnings during the comparison with samplesheet-provided values.

### `Added`

- The `seqkit sort` step prior to `bwa mem` alignment is now optional and disabled by default. This step was initially introduced to improve pipeline reproducibility. `bwa mem` is only deterministic when the order of the reads in its input fastq files remains unchanged, but since the order of the fastq files produced during `deacon` host removal depends on the relative speed of the different threads, this was not always the case.

### `Fixed`

- The resource requirements for `seqkit sort` were changed in order to (hopefully) avoid OOM (out-of-memory) errors when processing large fastq files when using the SLURM job scheduler. The process now requests more initial RAM and does not request additional CPU on retries (always 4 threads as per [the seqkit docs recommendations](https://bioinf.shenwei.me/seqkit/usage/#parallelization-of-cpu-intensive-jobs)).
- Fix bug in samplesheet parsing for single-end input where auto-detected flowcell information was not being assigned correctly to the meta map.
- Flowcell detection should now be compatible with (most) SRA/ENI fastq files (fastq header prefixed with space separated accession before the Illumina header).


## v0.1.0 - [2026-06-06]

Initial release of ITMmalaria/plasmovar, created with the [nf-core](https://nf-co.re/) template.

### `Added`

Added initial implementation of pipeline containing the following components:

0. Preparation of reference genome file(s) and regions (e.g., various indexes, bed/interval_lists, dictionaries, etc.)
1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Read trimming and adapter removal ([`fastp`](https://github.com/opengene/fastp))
3. Assess read composition ([`FastQ Screen`](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/))
4. Host read removal / decontamination ([`deacon`](https://github.com/bede/deacon) / ~~[`bbsplit`](https://sourceforge.net/projects/bbmap/)~~)
5. Read alignment against reference genome ([`bwa mem`](https://github.com/lh3/bwa) and [`seqkit`](https://bioinf.shenwei.me/seqkit/))
6. BAM file processing:
   - Duplicate marking and merging runs ([`gatk4 picard`](https://github.com/broadinstitute/gatk))
   - Coordinate sorting ([`samtools sort`](https://samtools.github.io/))
7. Summarise alignment statistics ([`samtools stats/flagstat`](https://samtools.github.io/) and [`mosdepth`](https://github.com/brentp/mosdepth))
8. Joint variant calling ([`gatk`](https://github.com/broadinstitute/gatk))
9. Variant filtering (variant quality score recalibration or hard filtering [`gatk`](https://github.com/broadinstitute/gatk))
10. Variant annotation ([`snpEff](https://pcingola.github.io/SnpEff/))
11. Create summary report ([`MultiQC`](http://multiqc.info/))

### `Fixed`

### `Dependencies`

### `Deprecated`
