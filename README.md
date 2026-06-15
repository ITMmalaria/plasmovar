# ITMmalaria/plasmovar

<!-- [![GitHub Actions CI Status](https://github.com/ITMmalaria/plasmovar/actions/workflows/nf-test.yml/badge.svg)](https://github.com/ITMmalaria/plasmovar/actions/workflows/nf-test.yml) -->
<!-- [![GitHub Actions Linting Status](https://github.com/ITMmalaria/plasmovar/actions/workflows/linting.yml/badge.svg)](https://github.com/ITMmalaria/plasmovar/actions/workflows/linting.yml) -->
<!-- [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX) -->
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)
[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.0.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.0.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**ITMmalaria/plasmovar** is a Nextflow pipeline for variant calling (snp/indel) of short-read whole-genome sequencing data of human malaria parasites (_Plasmodium_).



<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/community/brand/workflow-schematics#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

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
10. Variant annotation ([`snpEff`](https://pcingola.github.io/SnpEff/))
11. Create summary report ([`MultiQC`](http://multiqc.info/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

A minimal samplesheet for the pipeline would be a simple csv containing a `sample` and `fastq_1` column (and `fastq_2` for paired-end sequencing data). Each row represents a fastq file (single-end) or a pair of fastq files (paired end) derived from a particular biological sample.


```csv title="minimal-samplesheet.csv"
sample,fastq_1,fastq_2
sample-1,sample-1_R1.fastq.gz,sample-1_R2.fastq.gz
sample-2,sample-2_R1.fastq.gz,sample-2_R2.fastq.gz
```

More detailed information can be founded in the [usage documentation](./docs/usage.md).

Additionally required input files are reference genome fasta files ([see below](#reference-genome-files)) as well as other files depending on the selected options. For more info, please have a look at the [parameter documentation](./docs/parameters.md).

<!-- TODO: describe additional input files: reference fasta's, suggested versions, deacon index + construction, recommended parameters per species for e.g. intervals -->

Now, you can run the pipeline using:

```bash
NXF_SYNTAX_PARSER=v1 nextflow run ITMmalaria/plasmovar \
   -profile <docker/singularity/.../institute> \
  --input samplesheet.cvs \
  --reference_fasta PlasmoDB-68_Pfalciparum3D7_Genome.fasta \
  --skip_fastqscreen \
  --skip_hostremoval \
  --skip_annotation \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

> [!WARNING]
> Due to [various](https://nf-co.re/docs/developing/migration-guides/strict-syntax) [changes](https://docs.seqera.io/nextflow/strict-syntax) in Nextflow (specifically [string parameter changes](https://nf-co.re/blog/2026/parameter-types)), the pipeline currently only supports the legacy parser (e.g. export `NXF_SYNTAX_PARSER=v1`) or nextflow versions `<26.04` (e.g. `export NXF_VER=25.10.4`).

### Reference genome files

To run the pipeline, at minimum a reference genome for the main species of _Plasmodium_ must be provided. These can be obtained from e.g., [PlasmoDB](https://plasmodb.org/plasmo/app/downloads).

```yml
reference_fasta: "/home/pmoris/itg/datasets/reference-genomes/Pfalciparum/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"

```

Various (optional) components of the pipeline rely on additional reference genomes, e.g. for host read removal using `deacon` and assessing the read composition using `FastQ Screen`. These fasta files can be supplied using the following options in a `params.yml` file (or directly on the command-line).

```yml
hostremoval_reference: "./data/ref/GRCh38.chr21.fa.gz"
fastqscreen_fastas: "./data/ref/Pfalciparum/PlasmoDB-release-68/PlasmoDB-68_Pfalciparum3D7_Genome.fasta,./data/ref/GRCh38.chr21.fa.gz" # comma-separated list of all genomes to screen against
```

To use `gatk` base quality score recalibration and variant quality score recalibration, high-quality VCF files with known variants need to be supplied. For the former, a (subset) of variants passing quality filters in one of the MalariaGEN Pf# releases can be used ([https://www.malariagen.net/resource/34/](https://www.malariagen.net/resource/34/)), for the latter _in-vivo P. falciparum_ crosses between laboratory strains can be used ([available here](https://www.malariagen.net/data_package/pf-crosses-1-0/)), as described in the supplementary materials of the publication [Pf8: an open dataset of Plasmodium falciparum genome variation in 33,325 worldwide samples](https://doi.org/10.6084/m9.figshare.29153447).

For species without a well-described _truth_ set, a bootstrapping approach can be employed for BQSR instead (as described [here](https://gatk.broadinstitute.org/hc/en-us/articles/360035890831-Known-variants-Training-resources-Truth-sets)), although the need for BQSR is somewhat debated. For VQSR, the simplest alternative is to disable it and rely on hard filtering instead.

```yml
run_bqsr: true
bqsr_known_sites_vcf: "./data/3d7_hb3.combined.final.vcf.gz"
bqsr_known_sites_tbi: "./data/3d7_hb3.combined.final.vcf.gz.tbi"

vcf_filter_mode: 'vqsr'
vqsr_snp_resource_vcfs: [
    './data/7g8_gb4.combined.final.vcf.gz',
    './data/hb3_dd2.combined.final.vcf.gz',
    './data/3d7_hb3.combined.final.vcf.gz',
    ]
vqsr_snp_resource_labels: [
  '--resource:7g8_gb4,known=false,training=true,truth=true,prior=15.0',
  '--resource:hb3_dd2,known=false,training=true,truth=true,prior=15.0',
  '--resource:3d7_hb3,known=false,training=true,truth=true,prior=15.0',
  ]

vqsr_indel_resource_vcfs: [
    './data/7g8_gb4.combined.final.vcf.gz',
    './data/hb3_dd2.combined.final.vcf.gz',
    './data/3d7_hb3.combined.final.vcf.gz',
    ]
vqsr_indel_resource_labels: [
  '--resource:7g8_gb4,known=false,training=true,truth=true,prior=15.0',
  '--resource:hb3_dd2,known=false,training=true,truth=true,prior=15.0',
  '--resource:3d7_hb3,known=false,training=true,truth=true,prior=15.0',
]
```

> [!NOTE]
> In the future, we hope to further automate the pipeline so that relevant data files and recommended parameters (e.g. intervals) will be downloaded automatically by simply specifying the species of your samples.

## Credits

ITMmalaria/plasmovar was originally written by Pieter Moris at the Malariology Unit of the Institute of Tropical Medicine Antwerp.

We thank the Nextflow and nf-core community for their extensive assistance in the development of this pipeline, through documentation, guidelines, examples and the support provided via channels such as Slack and GitHub.

Parts of this pipeline were developed using suggestions from Seqera AI, but no code was copied verbatim without thorough manual curation and validation. Examples include e.g. best practices for channel manipulation, comparison with existing nf-core implementations validation, etc.

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
If you use ITMmalaria/plasmovar for your analysis, please cite it using the following doi: [pending].

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
