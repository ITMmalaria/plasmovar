# ITMmalaria/plasmovar: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.1.0 - [date]

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
