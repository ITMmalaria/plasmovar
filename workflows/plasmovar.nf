/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Quality control and pre-processing
include { addReadgroupToMeta                     } from '../subworkflows/local/utils_nfcore_plasmovar_pipeline'
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_TRIMMED               } from '../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_DECONTAMINATED        } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { FASTP                                  } from '../modules/nf-core/fastp/main'
include { FASTQSCREEN_BUILDFROMINDEX             } from '../modules/local/fastqscreen/buildfromindex/main'
include { FASTQSCREEN_FASTQSCREEN                } from '../modules/local/fastqscreen/fastqscreen/main'
// Host read removal
include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_INDEXER } from '../modules/nf-core/bbmap/bbsplit/main'
include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_MAPPER  } from '../modules/nf-core/bbmap/bbsplit/main'
include { DEACON_INDEX                           } from '../modules/nf-core/deacon/index/main'
include { DEACON_INDEX_DIFF                      } from '../modules/local/deacon/diff/main'
include { DEACON_FILTER                          } from '../modules/nf-core/deacon/filter/main'
// Prepare reference genome
include { CREATE_INTERVALS_BED                   } from '../modules/local/create_intervals_bed/main'
include { SAMTOOLS_FAIDX                         } from '../modules/nf-core/samtools/faidx/main'
include { BWA_INDEX                              } from '../modules/nf-core/bwa/index/main'
include { BWA_INDEX as BWA_INDEX_FASTQSCREEN     } from '../modules/nf-core/bwa/index/main'
// Alignment
include { SEQKIT_FASTQ_SORT                      } from '../subworkflows/local/seqkit_fastq_sort/main'
include { BWA_MEM                                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MARKDUP } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_STATS                         } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT                      } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS                      } from '../modules/nf-core/samtools/idxstats/main'
include { GATK4_MARKDUPLICATES                   } from '../modules/nf-core/gatk4/markduplicates/main'
include { MOSDEPTH                               } from '../modules/nf-core/mosdepth/main'
// Variant calling
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST                } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_INTERVALLISTTOOLS                } from '../modules/nf-core/gatk4/intervallisttools/main'
include { BQSR                                   } from '../subworkflows/local/bqsr/main.nf'
include { GATK4_HAPLOTYPECALLER                  } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT                 } from '../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                    } from '../modules/nf-core/gatk4/genotypegvcfs/main'
// Variant filtering
include { VARIANT_FILTERING_HARD                 } from '../subworkflows/local/variant_filtering_hard/main.nf'
include { VARIANT_FILTERING_VQSR                 } from '../subworkflows/local/variant_filtering_vqsr/main.nf'
// Variant annotation
include { SNPEFF_BUILD                           } from '../modules/local/snpeff/build/main'
include { SNPEFF_ANNOTATE                        } from '../modules/local/snpeff/annotate/main'
include { HTSLIB_BGZIPTABIX as BGZIP_SNPEFF_VCF  } from '../modules/nf-core/htslib/bgziptabix/main'
 // nf-core utils
include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_plasmovar_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PLASMOVAR {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()


    //
    // Parse skip_* and only_* options to decide which optional steps will be run
    //

    def only_flags = [
        only_index_reference : params.only_index_reference,
        only_hostremoval     : params.only_hostremoval,
        only_fastqscreen     : params.only_fastqscreen,
    ]
    def enabled_only_flags = only_flags.findAll { _k, v -> v }.keySet()
    if (enabled_only_flags.size() > 1) {
        log.error "The following --only_* options were used together: ${enabled_only_flags.join(', ')}"
        error "Stopping pipeline. Please select only one --only_* option."
    } else if (params.only_index_reference) {
        if (params.reference_index) {
            log.error "--only_index_reference option was selected, but a pre-built index was already supplied via --reference_index."
            error "Stopping pipeline. Please re-try with different options."
        }
        log.warn("--only_index_reference option was selected. All selected --skip_* options will be ignored and all optional steps will be skipped (qc, trimming, hostremoval, alignment, variant calling, annotation).")
        skip_qc             = true
        skip_trimming       = true
        skip_fastqscreen    = true
        skip_hostremoval    = true
        skip_alignment      = true
        skip_variantcalling = true
        skip_annotation     = true
    } else if (params.only_hostremoval) {
        log.warn("--only_hostremoval option was selected. All selected --skip_* options will be ignored and only quality control and host decontamination steps will be run (including reference indexing if required).")
        skip_qc             = false
        skip_trimming       = true
        skip_fastqscreen    = true
        skip_hostremoval    = false
        skip_alignment      = true
        skip_variantcalling = true
        skip_annotation     = true
    } else if (params.only_fastqscreen) {
        log.warn("--only_fastqscreen option was selected. All selected --skip_* options will be ignored and only FastQ Screen step will be run (including reference indexing if required).")
        skip_qc             = true
        skip_trimming       = true
        skip_fastqscreen    = false
        skip_hostremoval    = true
        skip_alignment      = true
        skip_variantcalling = true
        skip_annotation     = true
    } else {
        skip_qc             = params.skip_qc
        skip_trimming       = params.skip_trimming
        skip_fastqscreen    = params.skip_fastqscreen
        skip_hostremoval    = params.skip_hostremoval
        skip_alignment      = params.skip_alignment
        skip_variantcalling = params.skip_variantcalling
        skip_annotation     = params.skip_annotation
        if (skip_alignment && (!skip_variantcalling || !skip_annotation)) {
            log.warn("Since the --skip_alignment option was selected, the variant calling and annotation steps will also be skipped (even though either of --skip_variantcalling or --skip_annotation was not explicitly set).")
            skip_variantcalling = true
            skip_annotation     = true
        }
        if (skip_variantcalling && !skip_annotation) {
            log.warn("Since the --skip_variantcalling option was selected, the annotation step will also be skipped (even though --skip_annotation was not explictly set).")
            skip_annotation     = true
        }
    }

    // TODO: Allow index only option to build indices for all required steps (bwa alignment, host decontamination, etc.)
    // TODO: move reference genome + index building logic to separate subworkflows. Decide if everything should be bundled (easier for the only_index option) or included in the relevant step (e.g. build index for deacon in deacon section). Check how nf-core/rnaseq handles genome prep logic.

    // TODO: allow either a single fasta reference file to be supplied or multiple ones
    // via samplesheet column or via extra file listing species name + reference path
    // optionally providing pre-built indices (bwa, samtools, gatk, etc.)
    // for reference sheet example: https://nf-co.re/eager/dev/docs/usage#reference-input
    // also requires changes to bbsplit (needs both host and main reference)

    // TODO: test data -> see separate branch

    // TODO validate reference genomes + automatic download
    // // Validate genome selection
    // if (!params.genome) {
    //     error "Please specify a genome with --genome. Available: ${params.genomes.keySet().join(', ')}"
    // }
    // if (!params.genomes.containsKey(params.genome)) {
    //     error "Genome '${params.genome}' not found. Available: ${params.genomes.keySet().join(', ')}"
    // }

    //
    // Import input data from samplesheet
    //

    // Add read groups to meta
    if (!params.only_index_reference) {
        ch_fastq = ch_samplesheet.map { meta, fastqs -> addReadgroupToMeta(meta, fastqs) }
    }

    //
    // Prepare reference genome and associated files
    //

    // TODO: add to subworkflow? Add indexing for which downstream steps?

    // Create channel containing the reference fasta
    def ref = file(params.reference_fasta, checkIfExists: true)
    def ref_basename = ref.simpleName
    ch_ref_fasta = channel.value(
        [
            [id: ref_basename],
            ref,
        ]
    )

    // Index reference fasta as .fai
    SAMTOOLS_FAIDX(
        ch_ref_fasta.map{ meta, fasta -> [ meta, fasta, [] ]},
        []
    )
    ch_ref_fai = SAMTOOLS_FAIDX.out.fai.collect() // [[id:Pf3D7_01_v3], /path/to/ref.fa.gz.fai] - will be a value channel because its input is one too

    // create combined channel as input for various downstream processes
    // explicitly turn it into a value channel to avoid downstream issues in e.g. samtools_sort, even though println ch_ref_fasta.getClass() suggests it already is (groovyx.gpars.dataflow.DataflowVariable, instead of DataflowBroadcast)
    ch_ref_fasta_fai = ch_ref_fasta.join(ch_ref_fai).collect()    // [[id:Pf3D7_01_v3], /path/to/ref.fa.gz, /path/to/ref.fa.gz.fai]

    // Create or read bed file
    if (params.reference_bed) {
        bed_file = file(params.reference_bed, checkIfExists: true)
        ch_ref_bed = channel.value([[id: bed_file.simpleName], bed_file])
    } else {
        CREATE_INTERVALS_BED(ch_ref_fai)
        ch_versions = ch_versions.mix(CREATE_INTERVALS_BED.out.versions)
        ch_ref_bed = CREATE_INTERVALS_BED.out.bed   // [[id:PlasmoDB-68_Pfalciparum3D7_Genome], /path/to/work/e6/dd568ffd9b39576d3677b1518aa507/PlasmoDB-68_Pfalciparum3D7_Genome.bed]

        // TODO: prepare intervals alternative method to bin regions
        // https://nf-co.re/modules/gatk4_preprocessintervals/
        // https://gatk.broadinstitute.org/hc/en-us/articles/13832754597915-PreprocessIntervals
        // e.g.  generate consecutive bins of 1000 bases from the reference, useful for species with too many small regions in fasta
        // TODO add option to supply custom bed file (e.g. ampliseq)
        // TODO: sort bed file? https://bedtools.readthedocs.io/en/latest/content/tools/sort.html unix sort should be faster. Where is sorting important? Is it the order within regions or all regions as they appear in the fasta?
    }

    //
    // Perform initial FASTQC quality control on raw reads
    //

    if (!skip_qc) {
        FASTQC (
            ch_fastq
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.map{ _meta, file -> file })  // MultiQC's fastqc module requires the zip output - https://multiqc.info/modules/fastqc/
    }

    //
    // Trim reads using fastp
    //

    // TODO: optionally switch to subworkflow:
        // https://nf-co.re/subworkflows/fastq_trim_fastp_fastqc => could work as is
        // https://nf-co.re/subworkflows/fastq_fastqc_umitools_fastp => convert to local subworkflow and omit redundant UMI module?
        // allows to check for and filter out empty fastq files after trimming (or below threshold)
        // second version also collects list of empty fastq files and number of reads for multiqc report

        // https://github.com/OpenGene/fastp?tab=readme-ov-file#adapters
        // fastp module has --detect_adapter_for_pe option enabled, meaning it will first search for
        // adapters through overlap analysis on a per-read basis.
        // Auto-detection (based on first ~1M reads) is disabled by default for PE, but is re-enabled
        // in the fastp module via --detect_adapter_for_pe (as recommended by manual:
        //  For PE data, fastp will run a little slower if you specify the sequence adapters or enable
        //  adapter auto-detection, but usually result in a slightly cleaner output, since the overlap
        //  analysis may fail due to sequencing errors or adapter dimers.)
        // Order for PE reads: 1) overlap per-read, 2) auto-detection (based on subset),
        //                      3) --adapter_sequence_r1/2, 4) fasta file
        // ALSO CHANGE MODULES.CONFIG e.g.             withName: '.*:FASTQ_FASTQC_UMITOOLS_FASTP:FASTP'
    // TODO: add option for single-ended reads
    // TODO: check fastp on split fastq option for speed-up, as used by sarek: https://nf-co.re/sarek/3.4.2/docs/usage/#split-fastq-files
    if (params.trimming_before_hostremoval && !skip_trimming) {
        // create expected input for fastp module using read adapter file path from input parameters
        // channel: [ val(meta), [ path(reads) ], path(adapters) ]
        ch_fastp_input = ch_fastq
            .map { meta, reads -> [ meta, reads, params.fastp_adapter_fasta ?: [] ] }

        FASTP (
            ch_fastp_input,
            false,  // discard_trimmed_pass - Specify true to not write any reads that pass trimming thresholds. This can be used to use fastp for the output report only. Previously set to `!params.fastp_save_trimmed`
            params.fastp_save_trimmed_fail,
            params.fastp_save_merged
        )
        ch_fastq_trimmed = FASTP.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})  // MultiQC's fastp module relies only on the json output - https://multiqc.info/modules/fastp/

        // Re-run fastQC on trimmed reads
        if (!skip_qc) {
            FASTQC_TRIMMED(ch_fastq_trimmed)
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]})
        }

        // Rename output channel using the original generic name to pass it on to downstream processes, since this step was optional
        ch_fastq = ch_fastq_trimmed
    }

    //
    // Screen reads against different references using FastQ Screen to assess read composition
    //

    if (!skip_fastqscreen) {
        // Check if required options are provided
        if (!params.fastqscreen_index_dir && !params.fastqscreen_fastas) {
            log.error("Neither --fastqscreen_index_dir nor --fastqscreen_fastas were provided, but FastQ Screen step was selected.")
            error "Stopping pipeline. Please provide missing options."
        }

        // Create bwa indexes for each of the provided reference fastas
        if (!params.fastqscreen_index_dir) {

            // Parse comma-separated reference fasta files as input for FastQ Screen
            def ref_fastas = params.fastqscreen_fastas
                .split(',')
                .collect {
                    file(it.trim(), checkIfExists: true)
                }
            ch_fastqscreen_ref_fastas = Channel.fromList(
                ref_fastas.collect { fastqscreen_ref ->
                    [[id: fastqscreen_ref.baseName], fastqscreen_ref]
                }
            )
            // !NOTE: when using file.simpleName, all file suffixes are removed. In the case of a reference like GRCh38.chr21.fa.gz, this means that meta will hold GRCh38.bed_file. When this is passed through BWA_INDEX, the indices will be prefixed with file.baseName, i.e. GRCh38.chr21.fa.{amb,ann,bwt,pac,sa}, which causes a mismatch between the two.

            // Construct bwa indexes for the fastqscreen reference fastas, if they are not provided
            BWA_INDEX_FASTQSCREEN(ch_fastqscreen_ref_fastas)

            // Collect index directories into flattened list (genome names will be extracted from index filenames)
            ch_fastqscreen_indexes = BWA_INDEX_FASTQSCREEN.out.index.collect { _meta, index_dir -> index_dir }
        }
        // Use pre-supplied directory of indexes otherwise
        else {
            // Check for redundant input options
            if (params.fastqscreen_fastas) {
                log.warn("Both a prebuilt --fastqscreen_index_dir and a list of --fastqscreen_fastas fastas were provided; proceeding with pre-built index for FastQ Screen.")
            } else {
                log.info("Using pre-generated directory of indexes for FastQ Screen.")
            }

            fastqscreen_index_dir = file(params.fastqscreen_index_dir, checkIfExists: true)
            ch_fastqscreen_indexes = channel.fromPath("$fastqscreen_index_dir/*", type: "dir").collect()
        }

        // Create FastQ Screen configuration file and directory with index files
        FASTQSCREEN_BUILDFROMINDEX(ch_fastqscreen_indexes, "bwa")
        ch_versions = ch_versions.mix(FASTQSCREEN_BUILDFROMINDEX.out.versions.collect())
        database_ch = FASTQSCREEN_BUILDFROMINDEX.out.database

        // Run FastQ Screen on all input reads with the same database
        FASTQSCREEN_FASTQSCREEN(
            ch_fastq,
            database_ch,
            "bwa"
        )
        ch_versions = ch_versions.mix(FASTQSCREEN_FASTQSCREEN.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQSCREEN_FASTQSCREEN.out.txt.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(FASTQSCREEN_FASTQSCREEN.out.html.collect{it[1]})
    }

    //
    // Host read removal / decontamination
    //

    // TODO: add BWA decontamination option
    // TODO: move to subworkflow?
    // TODO: use file with list of fasta reference paths and names? Or samplesheet could mention which species to map against (in bbsplit (human) and fastq screen)
    // Examples:
    // https://github.com/nf-core/eager/blob/dev/modules/local/host_removal.nf
    // https://github.com/nf-core/taxprofiler/blob/1.1.7/subworkflows/local/shortread_hostremoval.nf

    if (!skip_hostremoval) {

        // use bbsplit for host read removal
        if (params.hostremoval_method == "bbsplit") {

            // Check if required options are provided
            if (!params.hostremoval_reference && !params.hostremoval_bbsplit_index) {
                log.error("Neither --hostremoval_bbsplit_index nor --hostremoval_reference were provided, but host removal with bbsplit was selected.")
                error "Stopping pipeline. Please provide missing options."
            }

            // create bbsplit index from host reference if not provided
            if (!params.hostremoval_bbsplit_index) {
                // Prepare channel with list of reference genomes to filter reads against.
                // Expected format for BBMAP_BBSPLIT module is:
                //      tuple val(other_ref_names), path (other_ref_paths)
                //      [['name'], [/path/to/fast.gz]]
                hostremoval_reference = file(params.hostremoval_reference, checkIfExists: true)
                ch_bbsplit_other_refs = channel
                    .value([[params.hostremoval_bbsplit_reference_name], [hostremoval_reference]])

                // create bbsplit index for filtering
                BBMAP_BBSPLIT_INDEXER (
                    [ [:], [] ],    // input reads can be omitted for building a new index
                    [],             // input index can be omitted for building a new index
                    channel.value(file(params.reference_fasta, checkIfExists: true)),   // primary reference
                    ch_bbsplit_other_refs,
                    true            // only perform index building step
                )
                ch_bbsplit_index = BBMAP_BBSPLIT_INDEXER.out.index
                // bbsplit.sh -Xmx6000M ref_primary="/path/to/primary_genome.fasta"  ref_human="/path/to/contaminant_genome.fa.gz" path=bbsplit_index_output threads=4
            }
            // Retrieve host index from input parameters otherwise
            else {
                // Check for redundant input options
                if (params.hostremoval_reference) {
                    log.warn("Both a prebuilt --hostremoval_bbsplit_index and a --hostremoval_reference fasta were provided; proceeding with pre-built index for bbsplit host read removal.")
                } else {
                    log.info("Using pre-built reference indexes for bbsplit host read removal.")
                }
                // Index needs to be the directory `bbsplit_index` which contains a ref subdir,
                // which in turn contains an index and genome subdir.
                ch_bbsplit_index = channel.value(file(params.hostremoval_bbsplit_index, checkIfExists: true))
            }

            // run bbsplit in map mode
            BBMAP_BBSPLIT_MAPPER (
                ch_fastq,           // input reads
                ch_bbsplit_index,   // index database
                [],                 // primary reference can be omitted if an index is provided
                [ [], [] ],         // other references can be omitted if an index is provided
                false               // do not build a new index
            )
            ch_fastq_hostremoved = BBMAP_BBSPLIT_MAPPER.out.primary_fastq
            ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBSPLIT_MAPPER.out.stats.collect{it[1]})
        }

        // use Deacon for host read removal
        // TODO: add prebuilt pangenome index https://github.com/bede/deacon?tab=readme-ov-file#prebuilt-indexes or refer to it in docs and host on zenodo
        // TODO: add deacon output to multiqc?
        else if (params.hostremoval_method == "deacon") {

            // Check if required options are provided
            if (!params.hostremoval_reference && !params.hostremoval_deacon_index) {
                log.error("Neither --hostremoval_deacon_index nor --hostremoval_reference were provided, but host removal with deacon was selected.")
                error "Stopping pipeline. Please provide missing options."
            }

            // Create deacon index from host reference if not provided
            if (!params.hostremoval_deacon_index) {
                def deacon_fasta = file(params.hostremoval_reference, checkIfExists: true)
                ch_deacon_fasta = channel.of(
                    tuple(
                        [ id: deacon_fasta.simpleName ],
                        file(deacon_fasta)
                    )
                )

                DEACON_INDEX(ch_deacon_fasta)
                ch_deacon_index = DEACON_INDEX.out.index
            }

            // Retrieve host index from input parameters otherwise
            else {
                // Check for redundant input options
                if (params.hostremoval_reference) {
                    log.warn("Both a prebuilt --hostremoval_deacon_index and a --hostremoval_reference fasta were provided; proceeding with pre-built index for deacon host read removal.")
                } else {
                    log.info("Using pre-built reference indexes for deacon host read removal.")
                }
                def deacon_index = file(params.hostremoval_deacon_index, checkIfExists: true)
                ch_deacon_index = channel.of(tuple([id: deacon_index.simpleName], file(deacon_index)))
            }

            // Subtract shared minimizers between parasite index and host index (see https://github.com/bede/deacon?tab=readme-ov-file#set-operations)
            if ( params.hostremoval_deacon_diff ) {
                log.info("Reference genome will be subtracted from pre-built reference index for deacon host removal.")
                DEACON_INDEX_DIFF(ch_ref_fasta, ch_deacon_index)
                ch_deacon_index = DEACON_INDEX_DIFF.out.index.map{ meta, index ->
                    def index_id = "${meta.id}_diff_${ref_basename}"
                    [ [ id: index_id ], index ]
                }
            }

            // Create input channel for deacon containing reads and index
            ch_deacon_input = ch_fastq
                .combine(ch_deacon_index)
                .map { meta, reads, _meta_index, index -> [ meta, index, reads] }

            // Filter reads using deacon against host index
            DEACON_FILTER(ch_deacon_input)
            ch_fastq_hostremoved = DEACON_FILTER.out.fastq_filtered
        }

        // Re-run fastQC on host-filtered reads
        if (!skip_qc) {
            FASTQC_DECONTAMINATED(ch_fastq_hostremoved)
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_DECONTAMINATED.out.zip.collect{it[1]})
        }
        // Rename output channel using the original generic name to pass it on to downstream processes, since this step was optional
        ch_fastq = ch_fastq_hostremoved
    }

    if (!params.trimming_before_hostremoval && !skip_trimming) {
        // create expected input for fastp module using read adapter file path from input parameters
        // channel: [ val(meta), [ path(reads) ], path(adapters) ]
        ch_fastp_input = ch_fastq
            .map { meta, reads -> [ meta, reads, params.fastp_adapter_fasta ?: [] ] }

        FASTP (
            ch_fastp_input,
            false,  // discard_trimmed_pass - Specify true to not write any reads that pass trimming thresholds. This can be used to use fastp for the output report only. Previously set to `!params.fastp_save_trimmed`
            params.fastp_save_trimmed_fail,
            params.fastp_save_merged
        )
        ch_fastq_trimmed = FASTP.out.reads
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})  // MultiQC's fastp module relies only on the json output - https://multiqc.info/modules/fastp/

        // Re-run fastQC on trimmed reads
        if (!skip_qc) {
            FASTQC_TRIMMED(ch_fastq_trimmed)
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMMED.out.zip.collect{it[1]})
        }

        // Rename output channel using the original generic name to pass it on to downstream processes, since this step was optional
        ch_fastq = ch_fastq_trimmed
    }

    //
    // Create bwa index for reference genome if index is not already provided
    //

    // TODO: move into subworkflow?

    // Construct bwa index for the reference fasta if it is not supplied by the user
    if (!params.reference_index) {
        BWA_INDEX (ch_ref_fasta)
        ch_bwa_index = BWA_INDEX.out.index.collect()    // collect() is required to create a re-usable value channel, otherwise it will only contain a single element which won't be emitted for each of the sample reads in ch_fastq
    } else {
        // If pre-made index is provided, check if it matches the supplied reference fasta
        // It should be a path to bwa directory containing *.{amb,ann,btw,pac,sa} files that share the same basename as the ref
        def index_dir = file(params.reference_index, checkIfExists: true)
        def required_extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']

        // Look for index files matching the reference basename
        def index_files = required_extensions.collectMany { ext ->
            file("${index_dir}/${ref_basename}*.${ext}", glob: true).toList()
        }

        // Warn if no index files are found at all
        if (index_files.isEmpty()) {
            error """
            Incomplete BWA index for ${ref} in ${index_dir}
            Expected: ${ref_basename}*.{amb,ann,bwt,pac,sa}
            """.stripIndent()
        }

        // Warn if there are multiple unique basenames
        // This is required because the bwa/mem module uses a find command to retrieve `amb` files
        def basenames = index_files.collect { file -> file.baseName }.unique()

        if (basenames.size() > 1) {
            error """
            Expected exactly one BWA index with basename '${ref_basename}' in ${index_dir}, but found: ${basenames.join(', ')}
            Matching files: ${index_files.collect { file -> "${file.name}" }.join(', ')}
            Please ensure only one set of index files matches the reference basename pattern.
            """.stripIndent()
        }

        def index_basename = basenames[0]

        // Check if all expected index files are present
        def missing_files = required_extensions.findAll { ext ->
            !file("${index_dir}/${index_basename}.${ext}").exists()
        }
        if (missing_files) {
            error """
            Incomplete BWA index for '${index_basename}' in ${index_dir}
            A complete BWA index requires all 5 files: ${index_basename}.{amb,ann,bwt,pac,sa}
            Missing files: ${missing_files.join(', ')}
            Please rebuild the index using: bwa index ${ref}
            """.stripIndent()
        }

        // Create the channel with validated index directory (value channel allows singleton to be reused)
        ch_bwa_index = channel.value([
            [ id: index_basename ],
            index_dir
        ])
    }

    //
    // bwa alignment to reference genome
    //

    if (!skip_alignment) {

        // Optionally sort fastq reads to make bwa deterministic - can have high memory requirements
        if (params.sort_fastq) {
            SEQKIT_FASTQ_SORT(ch_fastq)
            ch_fastq = SEQKIT_FASTQ_SORT.out.fastq_sorted
        }

        // Run bwa mem and sort resulting bam files
        sort_bam = true
        BWA_MEM (
            ch_fastq,
            ch_bwa_index,
            ch_ref_fasta,
            sort_bam        // markduplicates expects coordinate (or query) sorted input
        )
        ch_bam = BWA_MEM.out.bam

        // TODO: optionally enable CRAM output

        // TODO: avoid stalling when combining runs/lanes from the same sample/library
        // old sarek approach: https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L272
        // new sarek approach: https://github.com/nf-core/sarek/blob/20f41d1ce8b7ba296ee22adc71fe2da2ebcae93f/subworkflows/local/fastq_preprocess_gatk/main.nf#L163

        //
         // Mark duplicates using Picard and merge reads derived from same sample using RG values
        //
        // TODO: order bams to improve reproducibility for validation testing?

        // Group multi-lane/library samples by creating tuple with sample as grouping key
        ch_bam_grouped = ch_bam
            // remove unique distinguishing fields in meta data per sample
            .map { meta, bam ->
                [ meta - meta.subMap('id', 'read_group', 'single_end') + [ data_type: 'bam' ], bam ]
            }
            .groupTuple()
            .map { meta, bam -> [ [ id: meta.sample ] +  meta, bam ] }  // add new (sample-level) id

        // TODO make markdup optional for amplicon analysis => need to merge bams in another way
        // https://sites.google.com/a/broadinstitute.org/legacy-gatk-documentation/frequently-asked-questions/6057-At-what-point-should-I-merge-read-group-BAM-files-belonging-to-the-same-sample-into-a-single-file
        // MergeSamFiles (Picard)

        GATK4_MARKDUPLICATES(
            ch_bam_grouped,
            ch_ref_fasta.map{ _meta, fasta -> [ fasta ] }.collect(),
            ch_ref_fai.map{ _meta, fai -> fai }.collect()
        )
        ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.collect{it[1]})
        ch_bam_markdup = GATK4_MARKDUPLICATES.out.bam

        // TODO check sarek's method of collecting reports...or metrics? https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L447C1-L448C104

        // TODO: picard module could be an alternative, but it can only process 1 sample at a time
        // i.e. it does not merge samples while taking into account read groups.
        // Requires manual concatenation/combining (https://nf-co.re/modules/samtools_cat/ or
        // https://nf-co.re/modules/samtools_merge/), sorting
        // (https://nf-co.re/modules/samtools_sort/) and indexing
        // (https://nf-co.re/modules/samtools_index/) afterwards.
        // Possible advantage: multiple read groups linear scales ram requirements
        // https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR
        // PICARD_MARKDUPLICATES(
        //     ch_bam,
        //     ch_ref_fasta,
        //     ch_ref_fai
        // )

        //
        // samtools sort and index bam files and collect alignment stats
        //

        // TODO: move into subworkflow?
        // Examples:
        // https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_stats_samtools/main.nf
        // https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_markduplicates_samtools/main.nf
        // https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_sort_stats_samtools/main.nf
        // https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_stats_samtools/main.nf
        // https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/fastq_align_bwa/main.nf
        // https://github.com/nf-core/modules/blob/master/subworkflows/nf-core/bam_markduplicates_picard/main.nf

        // Sort and index duplicate marked bam files
        SAMTOOLS_SORT_MARKDUP(ch_bam_markdup, ch_ref_fasta_fai, 'bai')
        ch_bam_markdup_sort = SAMTOOLS_SORT_MARKDUP.out.bam
        ch_bam_bai = ch_bam_markdup_sort.join(SAMTOOLS_SORT_MARKDUP.out.index)

        SAMTOOLS_STATS(ch_bam_bai, ch_ref_fasta_fai)
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]})

        SAMTOOLS_FLAGSTAT(ch_bam_bai)
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{it[1]})

        SAMTOOLS_IDXSTATS(ch_bam_bai)
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS.out.idxstats.collect{it[1]})

        //
        // Module: mosdepth coverage statistics
        //

        ch_bam_bai_bed = ch_bam_bai.combine(ch_ref_bed.map { _meta, bed -> bed })
        MOSDEPTH(ch_bam_bai_bed, ch_ref_fasta, [])
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.per_base_bed.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_bed.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.quantized_bed.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.thresholds_bed.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]})

        // TODO: check behaviour of mosdepth when supplying a restricted list of intervals (e.g. ampliseq + QC metric for genes of interest during WGS)

        // TODO: check out other possible alignment statistics modules:
        // https://nf-co.re/modules/mosdepth/
        // https://nf-co.re/modules/samtools_depth/
        // https://nf-co.re/modules/samtools_coverage/
        // https://nf-co.re/modules/samtools_bedcov/
        // https://nf-co.re/modules/coverm_genome/
        // https://nf-co.re/modules/deeptools_bamcoverage/
        // https://nf-co.re/modules/bedtools_coverage/
        // https://nf-co.re/modules/deeptools_multibamsummary/
    }

    //
    // Call variants using GATK
    //

    if (!skip_variantcalling && !skip_alignment) {

        // TODO: refactor into subworkflow
        // TODO: compare with approach used here https://github.com/nf-core/genomicrelatedness/blob/dev/subworkflows/local/combine_cram_crai_intervals/main.nf

        //
        // Prepare intervals for scatter-gather processing
        //

        // flow:
        // BAM files (N samples × M intervals)
        //     ↓ GATK4_HAPLOTYPECALLER [parallel: N×M]
        // GVCF files (N samples × M intervals)
        //     ↓ Channel grouping
        // GVCF groups (M intervals, each with N samples)
        //     ↓ GATK4_GENOMICSDBIMPORT [parallel: M]
        // GenomicsDB (M databases, one per interval)
        //     ↓ GATK4_GENOTYPEGVCFS [parallel: M]
        // Variant filtering (M interval VCFs, each with all N samples)
        //     ↓ VARIANT_FILTERING subworkflow [parallel: M]
        // VCF files (M interval VCFs, each with all N samples)
        //     ↓ PICARD_MERGEVCFS
        // Final VCF (1 file with all N samples, all M intervals)

        // Process each sample x interval combination in parallel with HaplotpeCaller
        // Group GVCFs for all samples per interval for GenomicsDBImport
        // Process each interval in parallel with GenotypeGVCFs
        // Merge intervals

        // Create GATK dictionary
        GATK4_CREATESEQUENCEDICTIONARY(ch_ref_fasta)
        ch_ref_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

        // Convert BED to GATK IntervalList
        GATK4_BEDTOINTERVALLIST(ch_ref_bed, ch_ref_dict)

        // TODO: consider using https://gatk.broadinstitute.org/hc/en-us/articles/9570421542811-ScatterIntervalsByNs-Picard and https://nf-co.re/modules/picard_scatterintervalsbyns/

        // TODO: add simpler strategy that just splits bed file into 1 task per contig + add more meta data about contig name

        // Divide chrom/contigs (or bed targets) across interval_list files for scatter-gather parallel processing
        GATK4_INTERVALLISTTOOLS(GATK4_BEDTOINTERVALLIST.out.interval_list)
        ch_intervals = GATK4_INTERVALLISTTOOLS.out.interval_list
            // [ [interval_genome_meta.id], [interval_list, interval_list, ...] ]               (single element)
            .map { meta, interval_lists ->
                tuple(
                    meta,
                    interval_lists.withIndex()
                )
            }
            // [ [interval_genome_meta.id], [ [interval_list, 0], [interval_list, 1], ...] ]    (single element channel)
            .transpose()
            // [ [interval_genome_meta.id], [interval_list, 0] ]                                (multiple elements in channel, 1 for each interval_list)
            .map { meta, indexed_interval ->
                def (interval_list, idx) = indexed_interval
                def genome_id = meta.id // GATK4_INTERVALLISTTOOLS passes its input genome meta.id
                def interval_name = interval_list.simpleName
                tuple(
                    meta + [
                        id: "${genome_id}_${interval_name}",
                        genome_id: genome_id,
                        interval_index: idx,
                        interval_name: interval_name
                    ],
                    // meta + [
                    //     id: "${meta.id}_${interval_list.simpleName}",
                    //     genome_id: meta.id,
                    //     interval_index: idx,
                    //     interval_name: interval_list.simpleName
                    // ],
                    interval_list
                )
            }
            // [ [interval_genome_meta.id, interval_index, interval_name], interval_list ]      (multiple elements in channel, 1 for each interval_list)

        // Combine each sample's bam file with each interval_list
        // (or rather, copy the bam - bam files are not subset to the interval regions but passed in their entirety)
        ch_bam_bai_intervals = ch_bam_bai
            .combine(ch_intervals)
            .map { bam_meta, bam, bai, interval_meta, interval ->
                def combined_meta = bam_meta + [
                    id: "${bam_meta.id}_${interval_meta.interval_name}",    // add unique id per sample/interval combination
                    genome_id: interval_meta.genome_id,
                    interval_index: interval_meta.interval_index,
                    interval_name: interval_meta.interval_name
                ]
                [ combined_meta, bam, bai, interval ]
            }
            // [ [combined_meta.id, combined_meta.sample, combined_meta.num_entries, combined_meta.multiple_lanes, combined_meta.multiple_libraries, multiple_flowcells, combined_meta.data_type,
            //  combined_meta.genome_id, combined_meta.interval_index, combined_meta.interval_name],
            //  bam, bai, interval_list ]
            // (sample x interval)

        // TODO: check if bed or interval_list is preferred

        // Alternative approach to generate scattered interval
        // BIN_INTERVALS(ch_ref_bed, 1000000, [])
        // TODO sarek alternative approach: https://github.com/nf-core/sarek/blob/master/modules/local/create_intervals_bed/main.nf
        // https://github.com/nf-core/sarek/blob/master/subworkflows/local/prepare_intervals/main.nf
        // https://github.com/nf-core/sarek/blob/master/subworkflows/local/bam_variant_calling_haplotypecaller/main.nf

        // TODO: compare interval approach with https://github.com/nf-core/genomicrelatedness/blob/dev/modules/local/splitintervals/main.nf
        // seems to split multiple times inside each subworkflow

        // Base Quality Score Recalibration (BQSR)
        // TODO: should BQSR be enabled by default or not? See concerns in docs.

        if (params.run_bqsr ==  true) {
            // Check if required files are supplied
            if (!params.bqsr_known_sites_vcf || !params.bqsr_known_sites_tbi) {
                log.error("bqsr_known_sites_vcf and bqsr_known_sites_tbi were not provided, but they are required when enabling base quality score recalibration (bqsr).")
                error "Stopping pipeline. Please provide missing files."
            }

            if (params.bqsr_bed && !params.bqsr_scatter) {
                log.error("Custom BQSR intervals were provided, but scatter mode was disabled.")
                error "Stopping pipeline. Please enable scatter mode when providing BQSR intervals."
            }

            BQSR(
                ch_bam_bai_intervals,
                ch_ref_fasta,
                ch_ref_fai,
                ch_ref_dict,
                params.bqsr_known_sites_vcf,
                params.bqsr_known_sites_tbi,
                params.bqsr_scatter,
                params.bqsr_bed,
            )
            ch_multiqc_files = ch_multiqc_files.mix(BQSR.out.recalibration_table.collect{it[1]})

            // Prepare channel for HaplotypeCaller by scattering bam files over intervals again
            // because BQSR ApplyBQSR runs in non-scattered mode and drops interval-level metadata.
            ch_bam_bai_intervals = BQSR.out.bam_recalibrated
                .join(BQSR.out.bai_recalibrated, by: 0)
                .combine(ch_intervals)
                // re-add meta fields that were omitted during scatter-gather operations in BQSR
                .map { recalibrated_bam_meta, bam, bai, interval_meta, interval ->
                    def combined_meta = recalibrated_bam_meta + [
                        id: "${recalibrated_bam_meta.id}_${interval_meta.interval_name}",
                        genome_id: interval_meta.genome_id,
                        interval_index: interval_meta.interval_index,
                        interval_name: interval_meta.interval_name
                    ]
                    [ combined_meta, bam, bai, interval ]
                }

        // TODO: compare with two-pass approach https://github.com/nf-core/genomicrelatedness/blob/dev/subworkflows/local/base_quality_score_recalibration/main.nf
        // TODO: https://gatk.broadinstitute.org/hc/en-us/articles/360037433771-GatherBQSRReports
        // TODO: https://gatk.broadinstitute.org/hc/en-us/articles/360037066912-AnalyzeCovariates
        // cf. https://github.com/nf-core/genomicrelatedness/
        }

        // Call variants per sample per interval in parallel (sample x interval)
        GATK4_HAPLOTYPECALLER(
            ch_bam_bai_intervals.map{ combined_meta, bam, bai, interval ->
                [ combined_meta, bam, bai, interval, [] ] },
            ch_ref_fasta,
            ch_ref_fai,
            ch_ref_dict,
            [[:], []],  // dbsnp (optional)
            [[:], []]   // dbsnp_tbi (optional)
        )
        ch_gvcf = GATK4_HAPLOTYPECALLER.out.vcf      // [[meta], vcf.gz]
        ch_gvcf_tbi = GATK4_HAPLOTYPECALLER.out.tbi  // [[meta], tbi]

        // Gather all samples for each interval for genomicsDBimport
        ch_gvcf_by_interval = ch_gvcf                                       // [ [meta - sample + interval], gvcf ]                 (sample x interval)
            // Join gvcf and tbi by meta
            .join(ch_gvcf_tbi, by: 0)                                       // [ meta, gvcf, tbi ]                                  (sample x interval)
            .map { meta, gvcf, tbi ->
                def interval_key = "${meta.genome_id}_${meta.interval_name}"
                [ interval_key, [meta.sample, meta, gvcf, tbi] ]            // [ interval_key, [ meta.sample, meta, gvcf, tbi ] ]   (sample x interval)
            }
            // Group by interval
            .groupTuple(by: 0)                                              // [ interval_key, [ [ meta, gvcf, tbi ], ... ] ]       (1 per interval)
            .map { interval_key, entries ->
                def sorted = entries.sort { it[0] }
                def metas = sorted.collect { it[1] }
                def gvcfs = sorted.collect { it[2] }
                def tbis  = sorted.collect { it[3] }
                // Create new meta without sample info
                def interval_meta = [
                    id: interval_key,
                    genome_id: metas[0].genome_id,
                    interval_index: metas[0].interval_index,
                    interval_name: metas[0].interval_name
                ]
                [ interval_meta, gvcfs, tbis ]                              // [ interval_meta, gvcfs, tbis ]                       (1 per interval)
            }

        // Re-join inputs for GenomicsDB (= GVCFs for all samples grouped by interval) with the matching interval_list files
        ch_genomicsdb_input = ch_gvcf_by_interval                   // [ [interval_meta], [gvcfs], [tbis] ]                 (1 per interval)
            .map { interval_meta, gvcfs, tbis ->
                [ interval_meta.interval_name, interval_meta, gvcfs, tbis ]
            }
            .join(
                ch_intervals.map { interval_meta, interval_file ->  // [ interval_meta, interval_list ]                     (1 per interval)
                    [ interval_meta.interval_name, interval_file ]
                }
            )
            .map { _interval_name, meta, gvcfs, tbis, interval_file ->
                [ meta, gvcfs, tbis, interval_file, [], [] ]        // [ [interval_meta], [gvcfs], [tbis], interval_list ]  (sample x interval)
            }

        // Collect GVCFs for each sample into a GenomicsDB per interval
        GATK4_GENOMICSDBIMPORT(
            ch_genomicsdb_input,
            false,  // run_intlist
            false,  // run_updatewspace
            false   // input_map
        )
        ch_genomicsdb = GATK4_GENOMICSDBIMPORT.out.genomicsdb

        // Prepare channels for GenotypeGVCFs by matching GenomicsDB with corresponding interval_list file
        ch_genotype_input = ch_genomicsdb
            .map { meta, genomicsdb ->                          // [ [interval_meta], genomicsdb_dir ]                      (1 per interval)
                [ meta.interval_name, meta, genomicsdb ]
            }
            .join(
                ch_intervals.map { meta, interval_file ->       // [ [interval_meta], interval_list ]                       (1 per interval)
                    [ meta.interval_name, interval_file ]
                }
            )
            .map { _interval_name, meta, genomicsdb, interval_file ->
                [ meta, genomicsdb, [], interval_file, [] ]     // [ [interval_meta], genomicsdb, [], interval_list, [] ]   (sample x interval)
            }

        // Perform joint genotyping per interval
        GATK4_GENOTYPEGVCFS(
            ch_genotype_input,
            ch_ref_fasta,
            ch_ref_fai,
            ch_ref_dict,
            [[:], []],  // dbsnp (optional)
            [[:], []]   // dbsnp_tbi (optional)
        )
        ch_vcf_by_interval = GATK4_GENOTYPEGVCFS.out.vcf        // [[meta], vcf.gz]     (1 per interval)
        ch_vcf_tbi_by_interval = GATK4_GENOTYPEGVCFS.out.tbi    // [[meta], vcf.gz.tbi] (1 per interval)

        //
        // Variant filtering subworkflow (intervals)
        //

        if (params.vcf_filter_mode == 'hard') {
            //! TODO: check if order of intervals is as expected

            // TODO: optionally allow for fasta gzi index in case of bgzf compressed fasta file, see https://nf-co.re/modules/gatk4_variantfiltration. Needs to alter VARIANT_FILTERING_HARD input options to accept this
            // If your pipeline already has ch_gzi from reference preparation, use it.
            // // Otherwise, default to empty:
            // def ch_gzi = params.fasta.endsWith('.gz')
            //     ? channel.of([[id: 'genome'], file("${params.fasta}.gzi")])
            //     : [[:], []]

            VARIANT_FILTERING_HARD(
                ch_vcf_by_interval.join(ch_vcf_tbi_by_interval),    // [[meta], vcf, tbi]
                ch_ref_fasta,
                ch_ref_fai,
                ch_ref_dict,
                // ch_gzi,
            )
            ch_final_vcf = VARIANT_FILTERING_HARD.out.vcf_filter_added          // [[meta], gathered.vcf.gz]
            ch_final_vcf_tbi = VARIANT_FILTERING_HARD.out.vcf_filter_added_tbi  // [[meta], gathered.vcf.gz.tbi]

        } else if (params.vcf_filter_mode == 'vqsr' || params.vcf_filter_mode == 'VQSR') {
            // Some aspects of the VQSR approach were adapted from the MalariaGEN Pf8 pipeline: https://github.com/malariagen/malariagen-pf8-snp-indel-calling/blob/master/main.nf
            VARIANT_FILTERING_VQSR(
                ch_vcf_by_interval.join(ch_vcf_tbi_by_interval),    // [[meta], vcf, tbi]
                ch_ref_fasta,
                ch_ref_fai,
                ch_ref_dict
            )
            ch_final_vcf = VARIANT_FILTERING_VQSR.out.vcf_filter_added          // [[meta], gathered.vcf.gz]
            ch_final_vcf_tbi = VARIANT_FILTERING_VQSR.out.vcf_filter_added_tbi  // [[meta], gathered.vcf.gz.tbi]
        }
    }

    //
    // snpEff annotation
    //

    if (!skip_annotation && !skip_variantcalling && !skip_alignment) {
        // Validate required input options
        if (!params.reference_annotation) {
            log.error("No reference annotation file was provided, but snpEff option was selected.")
            error "Stopping pipeline. Please provide missing options."
        }

        ch_snpeff_input = ch_ref_fasta.map { meta_ref, fasta ->
            tuple(
                meta_ref,
                fasta,
                file(params.reference_annotation, checkIfExists: true),
                params.reference_cds ? file(params.reference_cds, checkIfExists: true) : [],
                params.reference_protein ? file(params.reference_protein, checkIfExists: true) : [],
            )
        }

        SNPEFF_BUILD(
            ch_snpeff_input,
            params.annotation_format ?: '',
            ch_ref_bed,
            channel.fromPath("$projectDir/assets/snpEff.config", checkIfExists: true),
            ch_ref_fasta.map { meta, _fasta -> meta.id },
            )

        SNPEFF_ANNOTATE(
            ch_final_vcf,
            SNPEFF_BUILD.out.db,
            SNPEFF_BUILD.out.config,
            ch_ref_fasta.map { meta, _fasta -> meta.id }
        )
        ch_multiqc_files = ch_multiqc_files.mix(SNPEFF_ANNOTATE.out.report.collect{it[1]})

        BGZIP_SNPEFF_VCF(
            SNPEFF_ANNOTATE.out.vcf.map { meta, vcf -> [meta, vcf, [], []] },
            "compress",
            true,
            "vcf"
        )
    }

    // TODO: add step for gatk VariantsToTable -V "${ann_dir}/${species}/combined.filter_added.ann.vcf" -F CHROM -F POS -F TYPE -GF GT -O "${ann_dir}/${species}/combined.filter_added.table"
    // TODO: add bcf normalizes as safety measure?
    // TODO: add CDS/region information to annotation as done by MalariaGEN Pf8? https://github.com/malariagen/malariagen-pf8-snp-indel-calling/blob/master/modules/variant_annotation.nf
    // TODO: check for other variant filtration options used in MalariaGEN Pf8? https://github.com/malariagen/malariagen-pf8-snp-indel-calling/blob/master/modules/variant_filtration.nf

    //
    // Collate and save software versions
    //

    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'plasmovar_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // Run MultiQC module
    //

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'plasmovar'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    emit:multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
