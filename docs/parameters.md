# ITMmalaria/plasmovar pipeline parameters

A short-read variant calling pipeline for malaria molecular surveillance

## Input/output options

Define where the pipeline should find input data and save output data.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `input` | Path to comma-separated file containing information about the samples in the experiment. <details><summary>Help</summary><small>You will need to create a design file with information about the samples in your experiment before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row.</small></details>| `string` |  | True |  |
| `outdir` | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` |  | True |  |
| `multiqc_title` | MultiQC report title. Printed as page header, used for filename if not otherwise specified. | `string` |  |  |  |
| `strict_lane_detection` | Handles behaviour when detecting and validating lane info for input fastq files. <details><summary>Help</summary><small>In strict mode, the pipeline will halt if 1) the lane provided in the samplesheet does not match the lane detected in the fastq file name, 2) the detected lanes do match between read pairs, and 3) if no lane was provided in the samplesheet and detection fails. Otherwise, only warnings are issued, except for read pair mismatches, which always cause an error. Automatic detection searches for file name patterns such as `_L001_`, `_L1_`, `.L01`, `_lane1_` or `_lane-01_`.</small></details>| `boolean` | True |  | True |
| `strict_flowcell_detection` | Whether to (attempt to) automatically extract flowcell IDs from from the input fastq files when this info is missing from the samplesheet. <details><summary>Help</summary><small>In strict mode, the pipeline will halt if 1) the flowcell provided in the samplesheet does not match the flowcell detected in the fastq header, 2) the detected flowcells do match between read pairs, and 3) if no flowcell was provided in the samplesheet and detection fails. Otherwise, only warnings are issued, except for read pair mismatches, which always cause an error. Attempts to parse the first FASTQ header and extract the flowcell ID. Expects CASAVA 1.8+ or the old Illumina format, and can also handle SRA/ENI accession prefix (space-separated).</small></details>| `boolean` | True |  | True |

## Reference genome options

Reference genome related files and options required for the workflow.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `reference_fasta` | Path to main reference genome in FASTA format. <details><summary>Help</summary><small>This parameter is *mandatory*. If you don't have a bwa index available for the main parasite genome, this will be generated for you automatically. Combine with `--save_index` to save bwa index for future runs.</small></details>| `string` |  | True |  |
| `reference_index` | Optional path to directory housing BWA index files. <details><summary>Help</summary><small>Expected files are `*.{amb,.ann,.bwt,.pac,.sa}`. Files should have the same name as the main parasite reference genome FASTA.</small></details>| `string` |  |  |  |
| `save_index` | Whether to save the generated bwa reference index. <details><summary>Help</summary><small>Set to true to save the generated bwa index for the main parasite reference. The bwa index can be supplied to future runs via `--reference_index` to avoid unnecessary recomputation.</small></details>| `boolean` | True |  |  |
| `reference_bed` | Optional path to a BED file containing the genomic regions to restrict variant calling to and to scatter/parallellize the main gatk variant calling steps by. <details><summary>Help</summary><small>If not supplied, all chromosomes/contigs (genomic regions) defined in the reference fasta will be used. These will be grouped into a manageable number of interval_lists to avoid problems with highly fragmented genomes with many short contigs. The number of which can be controlled via `scatter_count`. Note that a separate bed file can be provided for the BQSR step if needed.</small></details>| `string` |  |  |  |

## Skipping options and order of steps

Disable specific steps of the pipeline.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `skip_qc` | Skip QC steps. <details><summary>Help</summary><small>If enabled, all quality control steps will be skipped.</small></details>| `boolean` | False |  |  |
| `skip_trimming` | Skip fastq trimming steps. <details><summary>Help</summary><small>If enabled, read trimming step will be skipped.</small></details>| `boolean` | False |  |  |
| `skip_fastqscreen` | Skip FastQ Screen step. <details><summary>Help</summary><small>If enabled, reads will not be screened against different reference genomes to assess the composition.</small></details>| `boolean` | False |  |  |
| `skip_hostremoval` | Skip removal of human reads. <details><summary>Help</summary><small>If enabled, host reads will not be filtered out.</small></details>| `boolean` | False |  |  |
| `skip_alignment` | Skip alignment steps. <details><summary>Help</summary><small>If enabled, reads will not be aligned. Implies that any further downstream steps (variant calling and annotation) will be skipped as well.</small></details>| `boolean` | False |  |  |
| `skip_variantcalling` | Skip variant calling steps. <details><summary>Help</summary><small>If enabled, variants will not be called. Implies that any further downstream steps (annotation) will be skipped as well.</small></details>| `boolean` | False |  |  |
| `skip_annotation` | Skip variant annotation steps. <details><summary>Help</summary><small>If enabled, variants will not be annotated.</small></details>| `boolean` | False |  |  |
| `only_index_reference` | Only build reference indices and exit. <details><summary>Help</summary><small>If enabled, all optional steps will be skipped and only indices for reference genomes will be built.</small></details>| `boolean` | False |  |  |
| `only_hostremoval` | Only perform host removal steps. <details><summary>Help</summary><small>If enabled, only quality control and host decontamination will be performed (including building of reference index if required).</small></details>| `boolean` | False |  |  |
| `only_fastqscreen` | Only perform FastQ Screen steps. <details><summary>Help</summary><small>If enabled, only FastQ Screen step will be performed (including building of reference indces if required).</small></details>| `boolean` | False |  |  |
| `trimming_before_hostremoval` | Whether to perform trimming before or after host read removal. | `boolean` | True |  |  |

## Read trimming options for fastp

Configure fastq trimming options.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fastp_save_trimmed` | Save trimmed FastQ files. | `boolean` | False |  |  |
| `fastp_save_merged` | Save merged paired-end reads. <details><summary>Help</summary><small>See: https://github.com/OpenGene/fastp#merge-paired-end-reads</small></details>| `boolean` | False |  |  |
| `fastp_save_trimmed_fail` | Save reads that fail quality filters. <details><summary>Help</summary><small>See: https://github.com/OpenGene/fastp?tab=readme-ov-file#store-the-reads-that-fail-the-filters</small></details>| `boolean` | False |  |  |
| `fastp_adapter_fasta` | Optional path to FASTA file containing adapter sequences for trimming. <details><summary>Help</summary><small>fastp will automatically detect adapters if this is not supplied.</small></details>| `string` |  |  |  |
| `fastp_adapter_sequence` | Optional string to specify adapter for read 1. <details><summary>Help</summary><small>fastp will automatically detect adapters if this is not supplied.</small></details>| `string` |  |  |  |
| `fastp_adapter_sequence_r2` | Optional string to specify adapter for read 2. <details><summary>Help</summary><small>fastp will automatically detect adapters if this is not supplied.</small></details>| `string` |  |  |  |
| `fastp_overrepresentation` | Enable overrepresented sequence analysis (disabled by default). | `boolean` | False |  |  |
| `fastp_trim_poly_g` | Enable polyG tail trimming. <details><summary>Help</summary><small>When setting this to false, it will still be enabled by default for NextSeq/NovaSeq data. See: https://github.com/OpenGene/fastp?tab=readme-ov-file#polyg-tail-trimming. To completely disable polyG trimming for all inputs, pass `--disable_trim_poly_g` as an argument to `--fastp_extra_args`.</small></details>| `boolean` | True |  |  |
| `fastp_qualified_quality_phred` | The minimum quality value that a base needs to meet to be qualified for quality filtering of reads. Used in conjunction with --fastp_unqualified_percent_limit. Fastp default is >= Q15 when omitted. | `integer` |  |  |  |
| `fastp_unqualified_percent_limit` | Limits the percentage of bases that are allowed to be unqualified for quality filtering of reads (0-100). Fastp default is 40 percent when omitted. | `integer` |  |  |  |
| `fastp_length_required` | Reads shorter than this will be discarded. Fastp default is 15 bp when omitted. | `integer` |  |  |  |
| `fastp_cut_front` | Enable cutting by quality from the front (5' end) using a sliding window. | `boolean` | False |  |  |
| `fastp_cut_tail` | Enable cutting by quality from the tail (3' end) using a sliding window. | `boolean` | False |  |  |
| `fastp_cut_right` | Enable cutting by quality from the front (5' end) using a sliding window and dropping the window and everything to the right. | `boolean` | False |  |  |
| `fastp_disable_adapter_trimming` | Disable adapter trimming. By default this option is set to false, i.e. adapter trimming is performed. | `boolean` | False |  |  |
| `fastp_extra_args` | Extra arguments to pass to fastp. | `string` |  |  |  |

## FastQ Screen options

Configure options for fastq-screen.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `fastqscreen_fastas` | Comma separated list of file paths to reference genome fastas. | `string` |  |  |  |
| `fastqscreen_save_index` | Save FastQ Screen bwa index files and configuration file. | `boolean` | False |  |  |
| `fastqscreen_index_dir` | Optional path to directory with pre-generated bwa index files for FastQ Screen. | `string` |  |  |  |

## Host read removal options

Configure human host read removal options.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `hostremoval_save_filtered` | Save filtered (host-removed) FastQ files. | `boolean` | False |  |  |
| `hostremoval_save_index` | Save generated host removal index files. | `boolean` | False |  |  |
| `hostremoval_method` | Method to use for host read removal. (accepted: `deacon`\|`bbsplit`) <details><summary>Help</summary><small>Can be set to either `deacon` or `bbsplit`.</small></details>| `string` | deacon |  |  |
| `hostremoval_reference` | Path to host reference genome FASTA for read removal. | `string` |  |  |  |
| `hostremoval_bbsplit_index` | Path to BBSplit index directory for host removal. | `string` |  |  |  |
| `hostremoval_bbsplit_reference_name` | Reference name for BBSplit host removal. | `string` | human |  |  |
| `hostremoval_bbsplit_ambiguous2` | BBSplit ambiguous2 parameter - how to handle ambiguous reads. | `string` | all |  |  |
| `hostremoval_deacon_index` | Path to DEACON index directory for host removal. | `string` |  |  |  |
| `hostremoval_deacon_diff` | Whether to subtract the main reference genome index from the host filter index. | `boolean` | True |  |  |
| `hostremoval_deacon_abs_threshold` | Minimum absolute number of minimizer hits for a match in deacon filter. Deacon's default is 2 when omitted. | `integer` |  |  |  |

## Alignment options

Configure options for bwa alignment of (processed) reads against reference genome.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `save_bam_initial` | Whether to save the initial BAM files. <details><summary>Help</summary><small>These are the initial BAM files generated by BWA, prior to merging or duplicate marking.</small></details>| `boolean` | False |  |  |
| `save_bam` | Whether to save the final BAM files. <details><summary>Help</summary><small>These are the final BAM files generated by BWA, merged by sample (across lanes and libraries) and after optional duplicate marking. Note that this does not include any (optional) base quality score recalibration (see `save_bqsr_bam`).</small></details>| `boolean` | False |  |  |
| `sort_fastq` | Whether to sort fastq files before alignment. <details><summary>Help</summary><small>This can make the bam alignment files more deterministic, since the input order of reads in the fastq file has an effect and their output order is not guaranteed to be consistent during the read removal steps. However, the memory requirements can be very high, so this option is only recommended for testing purposes.</small></details>| `boolean` | False |  | True |

## GATK Read Group options

Set attributes for the GATK Read Groups.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `seq_center` | Sequencing center for GATK read group CN field. | `string` |  |  |  |
| `seq_platform` | Sequencing platform for GATK read group PL field. | `string` | ILLUMINA |  |  |

## GATK interval options

Configure intervals of genomic regions for scatter/gather parallelization of GATK.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `scatter_count` | Maximum number of interval_list files into which the provided intervals (provided by `reference_bed` or automatically extracted from the reference genome fasta) will be combined for the scattered/parallel processing of the GATK steps. | `integer` | 10 |  |  |
| `save_interval_scatter` | Save scattered interval lists (for debugging purposes). | `boolean` | False |  |  |

## GATK base quality score recalibration (BQSR) options

Configure options for BQSR.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `run_bqsr` | Whether to run BQSR or not. Disabled by default. | `boolean` | False |  |  |
| `bqsr_known_sites_vcf` | VCF file with known sites to be used for BQSR. | `string` |  |  |  |
| `bqsr_known_sites_tbi` | Index file for VCF with known sites to be used for BQSR. | `string` |  |  |  |
| `bqsr_scatter` | Whether to run BaseRecalibrator in scatter mode (= parallellized by intervals) or not. Enabled by default. <details><summary>Help</summary><small>Toggling this off (or providing alternative intervals for the BQSR step only via `bqsr_bed`) can avoid issues when the known sites VCF contains different genomic regions than those used by the other steps in the pipeline.</small></details>| `boolean` | True |  |  |
| `bqsr_bed` | Optional path to a BED file containing the intervals to scatter BQSR by and also restrict itself to. <details><summary>Help</summary><small>If not supplied, the intervals will be the same as those used by the other GATK steps (which are either based on the provided `reference_bed` or automatically extracted from the reference genome fasta).</small></details>| `string` |  |  |  |
| `save_bqsr_bam` | Whether to save recalibrated BAM files. Defaults to false. | `boolean` | False |  |  |

## Variant calling options

Configure GATK variant calling options

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `save_genomicsdb` | Save GenomicsDB workspace. | `boolean` | False |  |  |
| `save_gvcf` | Save intermediate GVCF files from HaplotypeCaller. Defaults to false. | `boolean` | False |  |  |
| `save_intermediate_vcf` | Save intermediate VCF files created during various variant calling and filtering steps. Defaults to false, which results in only the SNPEff annotated file being published. | `boolean` | False |  |  |
| `include_non_variant_sites` | Whether to emit every site that was jointly genotyped, including invariant sites. | `boolean` | True |  |  |
| `ploidy` | Sample ploidy to pass to GATK HaplotypeCaller. Defaults to GATK's 2 when omitted. | `integer` |  |  |  |

## Variant filtering options

Configure variant filtering tools.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `vcf_filter_mode` | VCF filtering mode. (accepted: `hard`\|`vqsr`\|`VQSR`) <details><summary>Help</summary><small>Can be set to either `hard` or `vqsr`.</small></details>| `string` | hard |  |  |
| `vcf_filter_snp_QD_filter` | QD (Quality by Depth) filter expression for SNPs. <details><summary>Help</summary><small>Filter expression for SNP QD metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | QD < 2.0 |  |  |
| `vcf_filter_snp_QD_name` | Filter name for QD filter. <details><summary>Help</summary><small>Name assigned to variants that fail the QD filter.</small></details>| `string` | QD-2 |  |  |
| `vcf_filter_snp_QUAL_filter` | QUAL (Quality) filter expression for SNPs. <details><summary>Help</summary><small>Filter expression for SNP QUAL metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` |  |  |  |
| `vcf_filter_snp_QUAL_name` | Filter name for QUAL filter. <details><summary>Help</summary><small>Name assigned to variants that fail the QUAL filter.</small></details>| `string` |  |  |  |
| `vcf_filter_snp_FS_filter` | FS (FisherStrand) filter expression for SNPs. <details><summary>Help</summary><small>Filter expression for SNP FS metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | FS > 60.0 |  |  |
| `vcf_filter_snp_FS_name` | Filter name for FS filter. <details><summary>Help</summary><small>Name assigned to variants that fail the FS filter.</small></details>| `string` | FS+60 |  |  |
| `vcf_filter_snp_ReadPosRankSum_filter` | ReadPosRankSum filter expression for SNPs. <details><summary>Help</summary><small>Filter expression for SNP ReadPosRankSum metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | ReadPosRankSum < -8.0 |  |  |
| `vcf_filter_snp_ReadPosRankSum_name` | Filter name for ReadPosRankSum filter. <details><summary>Help</summary><small>Name assigned to variants that fail the ReadPosRankSum filter.</small></details>| `string` | ReadPosRankSum-8 |  |  |
| `vcf_filter_snp_SOR_filter` | SOR (StrandOddsRatio) filter expression for SNPs. <details><summary>Help</summary><small>Filter expression for SNP SOR metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | SOR > 3.0 |  |  |
| `vcf_filter_snp_SOR_name` | Filter name for SOR filter. <details><summary>Help</summary><small>Name assigned to variants that fail the SOR filter.</small></details>| `string` | StrandOddsRatio+3 |  |  |
| `vcf_filter_snp_MQ_filter` | MQ (RMSMappingQuality) filter expression for SNPs. <details><summary>Help</summary><small>Filter expression for SNP MQ metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | MQ < 40.0 |  |  |
| `vcf_filter_snp_MQ_name` | Filter name for MQ filter. <details><summary>Help</summary><small>Name assigned to variants that fail the MQ filter.</small></details>| `string` | MQ-40 |  |  |
| `vcf_filter_snp_MQRankSum_filter` | MQRankSum filter expression for SNPs. <details><summary>Help</summary><small>Filter expression for SNP MQRankSum metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | MQRankSum < -12.5 |  |  |
| `vcf_filter_snp_MQRankSum_name` | Filter name for MQRankSum filter. <details><summary>Help</summary><small>Name assigned to variants that fail the MQRankSum filter.</small></details>| `string` | MQRankSum-12.5 |  |  |
| `vcf_filter_indel_QD_filter` | QD (Quality by Depth) filter expression for INDELs. <details><summary>Help</summary><small>Filter expression for INDEL QD metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | QD < 2.0 |  |  |
| `vcf_filter_indel_QD_name` | Filter name for INDEL QD filter. <details><summary>Help</summary><small>Name assigned to INDELs that fail the QD filter.</small></details>| `string` | QD-2 |  |  |
| `vcf_filter_indel_QUAL_filter` | QUAL (Quality) filter expression for INDELs. <details><summary>Help</summary><small>Filter expression for INDEL QUAL metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` |  |  |  |
| `vcf_filter_indel_QUAL_name` | Filter name for INDEL QUAL filter. <details><summary>Help</summary><small>Name assigned to INDELs that fail the QUAL filter.</small></details>| `string` |  |  |  |
| `vcf_filter_indel_FS_filter` | FS (FisherStrand) filter expression for INDELs. <details><summary>Help</summary><small>Filter expression for INDEL FS metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | FS > 200.0 |  |  |
| `vcf_filter_indel_FS_name` | Filter name for INDEL FS filter. <details><summary>Help</summary><small>Name assigned to INDELs that fail the FS filter.</small></details>| `string` | FS+200 |  |  |
| `vcf_filter_indel_ReadPosRankSum_filter` | ReadPosRankSum filter expression for INDELs. <details><summary>Help</summary><small>Filter expression for INDEL ReadPosRankSum metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | ReadPosRankSum < -20.0 |  |  |
| `vcf_filter_indel_ReadPosRankSum_name` | Filter name for INDEL ReadPosRankSum filter. <details><summary>Help</summary><small>Name assigned to INDELs that fail the ReadPosRankSum filter.</small></details>| `string` | ReadPosRankSum-20 |  |  |
| `vcf_filter_indel_SOR_filter` | SOR (StrandOddsRatio) filter expression for INDELs. <details><summary>Help</summary><small>Filter expression for INDEL SOR metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | SOR > 10.0 |  |  |
| `vcf_filter_indel_SOR_name` | Filter name for INDEL SOR filter. <details><summary>Help</summary><small>Name assigned to INDELs that fail the SOR filter.</small></details>| `string` | StrandOddsRatio+10 |  |  |
| `vcf_filter_indel_InbreedingCoeff_filter` | InbreedingCoefficient filter expression for INDELs. <details><summary>Help</summary><small>Filter expression for INDEL InbreedingCoeff metric during hard filtering. Set to empty string to disable this filter.</small></details>| `string` | InbreedingCoeff < -0.8 |  |  |
| `vcf_filter_indel_InbreedingCoeff_name` | Filter name for INDEL InbreedingCoeff filter. <details><summary>Help</summary><small>Name assigned to INDELs that fail the InbreedingCoeff filter.</small></details>| `string` | InbreedingCoeff-0.8 |  |  |
| `vqsr_snp_resource_vcfs` | List of VCF files to use as training resources for SNP VQSR. <details><summary>Help</summary><small>Array of paths to high-confidence VCF files used for training the SNP variant quality score recalibration model. Should be supplied in the same order as the resource labels. E.g. ['/absolute/path/to/7g8_gb4.combined.final.vcf.gz','/absolute/path/to/hb3_dd2.combined.final.vcf.gz','/absolute/path/to/3d7_hb3.combined.final.vcf.gz',]</small></details>| `array` |  |  |  |
| `vqsr_snp_resource_labels` | Resource labels for SNP VQSR training files. <details><summary>Help</summary><small>Array of GATK resource specifications for each SNP training VCF. Should be supplied in the same order as the resource VCFs. Format: '--resource:name,known=false,training=true,truth=true,prior=15.0'</small></details>| `array` |  |  |  |
| `vqsr_snp_annotations` | Annotations to use for SNP VQSR model training. <details><summary>Help</summary><small>Variant annotations to use for building the SNP VQSR model. To be supplied as copmlete GATK CLI arguments (e.g., '-an QD -an FS -an SOR').</small></details>| `string` | -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum |  |  |
| `vqsr_snp_tranches` | Tranches for SNP VQSR model. <details><summary>Help</summary><small>List of sensitivity tranches for the SNP VQSR model, in the format -tranche 100.0 -tranche 99.95 ... .</small></details>| `string` | -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 |  |  |
| `vqsr_indel_resource_vcfs` | List of VCF files to use as training resources for INDEL VQSR. <details><summary>Help</summary><small>Array of paths to high-confidence VCF files used for training the INDEL variant quality score recalibration model. Should be supplied in the same order as the resource labels.</small></details>| `array` |  |  |  |
| `vqsr_indel_resource_labels` | Resource labels for INDEL VQSR training files. <details><summary>Help</summary><small>Array of GATK resource specifications for each INDEL training VCF. Should be supplied in the same order as the resource VCFs. Format: '--resource:name,known=false,training=true,truth=true,prior=15.0'</small></details>| `array` |  |  |  |
| `vqsr_indel_annotations` | Annotations to use for INDEL VQSR model training. <details><summary>Help</summary><small>Variant annotations to use for building the INDEL VQSR model. To be supplied as complete GATK CLI arguments (e.g., '-an QD -an FS -an SOR').</small></details>| `string` | -an QD -an FS -an SOR -an MQRankSum -an ReadPosRankSum |  |  |
| `vqsr_indel_tranches` | Tranches for INDEL VQSR model. <details><summary>Help</summary><small>List of sensitivity tranches for the INDEL VQSR model, in the format -tranche 100.0 -tranche 99.95 ... .</small></details>| `string` | -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.8 -tranche 99.6 -tranche 99.5 -tranche 99.4 -tranche 99.3 -tranche 99.0 -tranche 98.0 -tranche 97.0 -tranche 90.0 |  |  |
| `vqsr_snp_max_gaussians` | Max number of Gaussians for the positive model when recalibrating SNPs. <details><summary>Help</summary><small>When not provided, the GATK default of 8 will be used.</small></details>| `integer` |  |  |  |
| `vqsr_indel_max_gaussians` | Max number of Gaussians for the positive model when recalibrating INDELs. <details><summary>Help</summary><small>When not provided, the GATK default of 8 will be used.</small></details>| `integer` |  |  |  |
| `vqsr_snp_truth_sensitivity_filter_level` | The truth sensitivity level at which to start filtering during GATK ApplyVQSR in SNP mode. <details><summary>Help</summary><small>Variants in tranches that fall below the specified truth sensitivity filter level have their FILTER field annotated with the corresponding tranche level.</small></details>| `number` | 99.5 |  |  |
| `vqsr_indel_truth_sensitivity_filter_level` | The truth sensitivity level at which to start filtering during GATK ApplyVQSR in INDEL mode. <details><summary>Help</summary><small>Variants in tranches that fall below the specified truth sensitivity filter level have their FILTER field annotated with the corresponding tranche level.</small></details>| `number` | 99.5 |  |  |

## Annotation options

Configure snpEff for annotation of variants.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `reference_annotation` | Path to gene annotation file (gff3/gtf) used for building the snpEff database. | `string` |  |  |  |
| `reference_cds` | Path to CDS FASTA file to use as a check during the snpEff database building step. | `string` |  |  |  |
| `reference_protein` | Path to protein FASTA file to use as a check during the snpEff database building step. | `string` |  |  |  |
| `annotation_format` | Type of annotation file (gff or gtf). (accepted: `gff`\|`gff3`\|`gtf`) | `string` |  |  |  |
| `save_snpeff_db` | Save snpEff database of the reference genome. | `boolean` | False |  |  |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `custom_config_version` | Git commit id for Institutional configs. | `string` | master |  | True |
| `custom_config_base` | Base directory for Institutional configs. <details><summary>Help</summary><small>If you're running offline, Nextflow will not be able to fetch the institutional config files from the internet. If you don't need them, then this is not a problem. If you do need them, you should download the files from the repo and tell Nextflow where to find them with this parameter.</small></details>| `string` | https://raw.githubusercontent.com/nf-core/configs/master |  | True |
| `config_profile_name` | Institutional config name. | `string` |  |  | True |
| `config_profile_description` | Institutional config description. | `string` |  |  | True |
| `config_profile_contact` | Institutional config contact information. | `string` |  |  | True |
| `config_profile_url` | Institutional config URL link. | `string` |  |  | True |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter | Description | Type | Default | Required | Hidden |
|-----------|-----------|-----------|-----------|-----------|-----------|
| `version` | Display version and exit. | `boolean` |  |  | True |
| `publish_dir_mode` | Method used to save pipeline results to output directory. (accepted: `symlink`\|`rellink`\|`link`\|`copy`\|`copyNoFollow`\|`move`) <details><summary>Help</summary><small>The Nextflow `publishDir` option specifies which intermediate files should be saved to the output directory. This option tells the pipeline what method should be used to move these files. See [Nextflow docs](https://www.nextflow.io/docs/latest/process.html#publishdir) for details.</small></details>| `string` | copy |  | True |
| `max_multiqc_email_size` | File size limit when attaching MultiQC reports to summary emails. | `string` | 25.MB |  | True |
| `monochrome_logs` | Do not use coloured log outputs. | `boolean` |  |  | True |
| `multiqc_config` | Custom config file to supply to MultiQC. | `string` |  |  | True |
| `multiqc_logo` | Custom logo file to supply to MultiQC. File name must also be set in the MultiQC config file | `string` |  |  | True |
| `multiqc_methods_description` | Custom MultiQC yaml file containing HTML including a methods description. | `string` |  |  |  |
| `validate_params` | Boolean whether to validate parameters against the schema at runtime | `boolean` | True |  | True |
| `pipelines_testdata_base_path` | Base URL or local path to location of pipeline test dataset files | `string` | https://raw.githubusercontent.com/nf-core/test-datasets/ |  | True |
| `trace_report_suffix` | Suffix to add to the trace report filename. Default is the date and time in the format yyyy-MM-dd_HH-mm-ss. | `string` |  |  | True |
| `help` | Display the help message. | `['boolean', 'string']` |  |  |  |
| `help_full` | Display the full detailed help message. | `boolean` |  |  |  |
| `show_hidden` | Display hidden parameters in the help message (only works when --help or --help_full are provided). | `boolean` |  |  |  |
