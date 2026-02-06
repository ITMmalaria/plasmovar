/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Quality control and pre-processing
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
include { BWA_MEM                                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX                         } from '../modules/nf-core/samtools/index/main'
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
include { BIN_INTERVALS                          } from '../modules/local/bin_intervals/main'
include { GATK4_HAPLOTYPECALLER                  } from '../modules/nf-core/gatk4/haplotypecaller/main'
include { GATK4_GENOMICSDBIMPORT                 } from '../modules/nf-core/gatk4/genomicsdbimport/main'
include { GATK4_GENOTYPEGVCFS                    } from '../modules/nf-core/gatk4/genotypegvcfs/main'
include { GATK4_MERGEVCFS                        } from '../modules/nf-core/gatk4/mergevcfs/main'
include { TABIX_TABIX                            } from '../modules/nf-core/tabix/tabix/main'
// Variant annotation
include { SNPEFF_BUILD                           } from '../modules/local/snpeff/build/main'
include { SNPEFF_ANNOTATE                        } from '../modules/local/snpeff/annotate/main'
include { paramsSummaryMap                       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_plasmovar_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE ADDITIONAL FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: move to samplesheet prep subworkflow?

// Add GATK read group to meta and remove lane
// Adapted from nf-core/sarek: https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L940
// See read group info here: https://gatk.broadinstitute.org/hc/en-us/articles/360035890671-Read-groups
// List of RG fields is defined by SAM format: https://samtools.github.io/hts-specs/SAMv1.pdf
// Recommended usage: https://support.sentieon.com/appnotes/read_groups/
def addReadgroupToMeta(meta, files) {
    // Note: needs to be run on initial file path before processing,
    // otherwise the path in the work dir will not be found
    def flowcell = flowcellLaneFromFastq(files[0])
    if ( !meta.single_end && flowcell != flowcellLaneFromFastq(files[1]) ){
        error("Flowcell ID does not match for paired reads of sample ${meta.id} - ${files}")
    }

    // Define RG values, only adding non-empty values
    // This avoids problems with Picard and the need for the `VALIDATION_STRINGENCY LENIENT` option,
    // since there will not be any empty fields like CN in the following example:
    // @RG	ID:106264-002-079.22NY35LT3.L007	CN:	PU:22NY35LT3.L007	SM:106264-002-079	LB:106264-002-079	PL:ILLUMINA

    // If we cannot read the flowcell ID from the fastq file, then we don't use it
    def RG_ID = flowcell ? "${meta.sample}.${flowcell}.${meta.lane}" : "${meta.sample}.${meta.lane}"

    // For LB, check if library is defined in input samplesheet (trim handles case of empty string)
    def RG_LB = meta.library?.trim() ? "${meta.sample}.${meta.library}" : "${meta.sample}"

    def RG_map = [
        ID: RG_ID,
        SM: meta.sample,
        LB: RG_LB,
        PL: params.seq_platform,                            // defaults to ILLUMINA (nextflow.config)
        CN: params.seq_center?.trim() ?: null,              // most likely blank
        PU: flowcell ? "${flowcell}.${meta.lane}" : null,   // used by BQSR if present, otherwise defaults to using ID
    ]

    // Build RG string: filter out empty fields and format as RG string `key:value`
    def RG_fields = RG_map
            // example input map: [ID: "sample1", SM: "sample1", CN: null, PU: "flowcell.L001"]
        .findAll { _key, value ->
            value != null && value != ''
        }   // returns map: [ID: "sample1", SM: "sample1", PU: "flowcell.L001"]
        .collect { key, value ->
            "${key}:${value}"
        }   // returns list: ["ID:sample1", "SM:sample1", "PU:flowcell.L001"]

    def read_group = "\"@RG\\t${RG_fields.join('\\t')}\""

    meta  = meta - meta.subMap('lane', 'library') + [read_group: read_group.toString()]

    return [ meta, files ]
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
// Adapted from nf-core/sarek https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L953
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// TODO: see rnaseq example bbsplit fasta list

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PLASMOVAR {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()

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

    // TODO create channel holding reference genome
    // Note: bwa_index expects meta channel, whereas bbsplit_indexer just needs a path
    // ch_reference = channel.value(file(params.reference_fasta))
    // ch_reference = channel.fromPath(params.reference_fasta).map{ref -> tuple (ref.simpleName, ref)}
    // reference = [[ id:'reference', primary:true ], file(params.reference_fasta)]

    // TODO: allow either a single fasta reference file to be supplied or multiple ones
    // via samplesheet column or via extra file listing species name + reference path
    // optionally providing pre-built indices (bwa, samtools, gatk, etc.)
    // for reference sheet example: https://nf-co.re/eager/dev/docs/usage#reference-input
    // also requires changes to bbsplit (needs both host and main reference)

    // TODO: integrate all modules/subworkflows into multiQC: https://nf-co.re/docs/contributing/tutorials/adding_modules_to_pipelines

    // TODO: https://nf-co.re/modules/fastqscreen_fastqscreen

    // TODO: describe samplesheet columns in readme as done here https://nf-co.re/sarek/3.4.4/docs/usage/

    // TODO: test data -> see separate branch

    // TODO validate reference genomes + automatic download
    // // Validate genome selection
    // if (!params.genome) {
    //     error "Please specify a genome with --genome. Available: ${params.genomes.keySet().join(', ')}"
    // }
    // if (!params.genomes.containsKey(params.genome)) {
    //     error "Genome '${params.genome}' not found. Available: ${params.genomes.keySet().join(', ')}"
    // }

    // // Create reference channels
    // ch_fasta = channel.value([
    //     [id: params.genome],
    //     file(params.fasta)
    // ])


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
    // ch_ref_fasta = channel.fromPath(params.reference_fasta)
    //     .map { ref -> tuple([ id: ref.simpleName ], ref) }
    //     .collect()   // or .first()
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
        ch_ref_fasta,
        [[:], []],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_fai      = SAMTOOLS_FAIDX.out.fai    // [[id:Pf3D7_01_v3], /path/to/ref.fa.gz.fai]

    // Create or read bed file
    if (params.reference_bed) {
        bed_file = file(params.reference_bed, checkIfExists: true)
        ch_ref_bed = channel.value([[id: bed_file.simpleName], bed_file])
    } else {
        CREATE_INTERVALS_BED(ch_fai)
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
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})  // MultiQC's fastqc module requires the zip output - https://multiqc.info/modules/fastqc/
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
    // TODO: check fastp on split fastq option: https://nf-co.re/sarek/3.4.2/docs/usage/#split-fastq-files
    if (!skip_trimming) {
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
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
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
    // Screen reads against different references using FastQ Screen to assess composition
    //

    if (!skip_fastqscreen) {
        // TODO: add validation that at least one of fastqscreen_index_dir or fastqscreen_fastas needs to be supplied + add precedence or make mutually exclusive.
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
            log.info("Using pre-generated directory of indexes for FastQ Screen.")
            fastqscreen_index_dir = file(params.fastqscreen_index_dir, checkIfExists: true)
            ch_fastqscreen_indexes = channel.fromPath("$fastqscreen_index_dir/*", type: "dir").collect()
        }

        // Create FastQ Screen configuration file and directory with index files
        FASTQSCREEN_BUILDFROMINDEX(ch_fastqscreen_indexes, "bwa")
        ch_versions = ch_versions.mix(FASTQSCREEN_BUILDFROMINDEX.out.versions.first())
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
    // TODO: use file with list of fasta reference paths and names?
    // samplesheet could contain species
    // cf. samplesheet could mention which species to map against in bbsplit (human) and fastq screen
    // TODO: examples:
    // https://github.com/nf-core/eager/blob/dev/modules/local/host_removal.nf
    // https://github.com/nf-core/taxprofiler/blob/1.1.7/subworkflows/local/shortread_hostremoval.nf

    if (!skip_hostremoval) {

        // use bbsplit for host read removal
        if (params.hostremoval_method == "bbsplit") {
            // use provided BBSplit index if supplied or generate from scratch otherwise
            if (!params.hostremoval_bbsplit_index) {
                // Prepare channel with list of reference genomes to filter reads against.
                // Expected format for BBMAP_BBSPLIT module is:
                //      tuple val(other_ref_names), path (other_ref_paths)
                //      [['name'], [/path/to/fast.gz]]
                channel.from( [
                    [params.hostremoval_bbsplit_reference_name,
                    params.hostremoval_reference]
                ] )
                    .collect{ id, fasta -> [ [id], [file(fasta, checkIfExists: true)] ] }
                    .set { ch_bbsplit_other_refs }

                // create bbsplit index for filtering
                BBMAP_BBSPLIT_INDEXER (
                    [ [:], [] ],
                    [],
                    channel.value(file(params.reference_fasta, checkIfExists: true)),
                    ch_bbsplit_other_refs,
                    true
                )
                ch_bbsplit_index = BBMAP_BBSPLIT_INDEXER.out.index
                ch_versions = ch_versions.mix(BBMAP_BBSPLIT_INDEXER.out.versions.first())
                // bbsplit.sh -Xmx6000M ref_primary="/path/to/primary_genome.fasta"  ref_human="/path/to/contaminant_genome.fa.gz" path=bbsplit_index_output threads=4
            } else {
                // Index needs to be the directory `genome/index/bbsplit` which contains a ref subdir,
                // which in turn contains an index and genome subdir.
                log.info("Using pre-supplied reference fasta for bbsplit host removal")
                ch_bbsplit_index = channel.value(file(params.hostremoval_bbsplit_index, checkIfExists: true))
            }

            // run bbsplit in map mode
            BBMAP_BBSPLIT_MAPPER (
                ch_fastq,
                ch_bbsplit_index,
                [],
                [ [], [] ],
                false
            )
            ch_fastq_hostremoved = BBMAP_BBSPLIT_MAPPER.out.primary_fastq
            ch_versions = ch_versions.mix(BBMAP_BBSPLIT_MAPPER.out.versions.first())
            // ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBSPLIT_MAPPER.out.stats.collect{it[1]})
            // TODO multiqc bbsplit not showing up due to bug https://github.com/MultiQC/MultiQC/pull/1513
        }

        // use Deacon for host read removal
        // TODO: add prebuilt pangenome index https://github.com/bede/deacon?tab=readme-ov-file#prebuilt-indexes
        // TODO: add deacon output to multiqc?
        else if (params.hostremoval_method == "deacon") {

            // create Deacon index if not provided
            if (!params.hostremoval_deacon_index) {
                def deacon_fasta = file(params.hostremoval_reference, checkIfExists: true)
                ch_deacon_fasta = channel.of(
                    tuple(
                        [ id: deacon_fasta.simpleName ],
                        file(deacon_fasta)
                    )
                )

                DEACON_INDEX(ch_deacon_fasta)
                ch_versions = ch_versions.mix(DEACON_INDEX.out.versions.first())
                ch_deacon_index = DEACON_INDEX.out.index

            } else {
                // retrieve index from input parameters
                def deacon_index = file(params.hostremoval_deacon_index, checkIfExists: true)
                ch_deacon_index = channel.of(
                    tuple(
                        [ id: deacon_index.simpleName ],
                        file(deacon_index)
                    )
                )
            }

            if ( params.hostremoval_deacon_diff ) {
                DEACON_INDEX_DIFF(ch_ref_fasta, ch_deacon_index)
                ch_versions = ch_versions.mix(DEACON_INDEX_DIFF.out.versions.first())
                ch_deacon_index = DEACON_INDEX_DIFF.out.index.map{ meta, index ->
                    def index_id = "${meta.id}_diff_${ref_basename}"
                    [ [ id: index_id ], index ]
                }
            }

            // Create input channel for deacon containing reads and index
            ch_deacon_input = ch_fastq
                .combine(ch_deacon_index)
                .map { meta, reads, _meta_index, index -> [ meta, index, reads] }

            // filter reads using deacon against host index
            DEACON_FILTER(ch_deacon_input)
            ch_fastq_hostremoved = DEACON_FILTER.out.fastq_filtered
            ch_versions = ch_versions.mix(DEACON_FILTER.out.versions.first())
        }

        // Re-run fastQC on host-filtered reads
        if (!skip_qc) {
            FASTQC_DECONTAMINATED(ch_fastq_hostremoved)
            ch_multiqc_files = ch_multiqc_files.mix(FASTQC_DECONTAMINATED.out.zip.collect{it[1]})
        }
        // Rename output channel using the original generic name to pass it on to downstream processes, since this step was optional
        ch_fastq = ch_fastq_hostremoved
    }

    //
    // Create bwa index for reference genome if index is not already provided
    //

    // TODO: move into subworkflow

    // Construct bwa index for the reference fasta if it is not supplied by the user
    if (!params.reference_index) {
        BWA_INDEX (ch_ref_fasta)
        ch_bwa_index = BWA_INDEX.out.index.collect()    // collect() is required to create a re-usable value channel, otherwise it will only contain a single element which won't be emitted for each of the sample reads in ch_fastq
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    } else {
        // If pre-made index is provided, check if it matches the supplied reference
        // It should be a path to bwa directory containing *.{amb,ann,btw,pac,sa} files

        def index_dir = file(params.reference_index, checkIfExists: true)
        def required_extensions = ['amb', 'ann', 'bwt', 'pac', 'sa']

        // Collect index files matching the reference basename
        // def index_files = required_extensions.collect { ext ->
        //     file("${index_dir}/${ref_basename}*.${ext}", glob: true)
        //     // glob is required to handle names like reference.fa.amb vs reference.amb
        // }

        def index_files = required_extensions.collectMany { ext ->
            def matches = file("${index_dir}/${ref_basename}*.${ext}", glob: true)
            // glob is required to handle names like reference.fa.amb vs reference.amb
            // in case of no glob results, we still need to return an (empty) list
            matches instanceof List ? matches : [matches]
        }

        // println "DEBUG: Found ${index_files.size()} files with glob pattern ${ref_basename}*"
            // index_files.each { println "DEBUG:   - ${it.name}" }

        if (index_files.isEmpty()) {
            error """
            No BWA index files found for reference '${ref}' in ${index_dir}

            Expected files matching pattern: ${ref_basename}*.{amb,ann,bwt,pac,sa}

            Please ensure:
            1. The index directory contains BWA index files
            2. Index file basenames start with: ${ref_basename}
            3. Index was built using: bwa index ${ref}
            """.stripIndent()
        }

        // extract index basename
        def basenames = index_files.collect { it.baseName }.unique()

        // println "DEBUG: Unique basenames found: ${basenames}"

        if (basenames.size() > 1) {
            error """
            Multiple index file basenames found matching '${ref_basename}*' in ${index_dir}
            Found basenames: ${basenames.join(', ')}

            Matching files:
            ${index_files.collect { "  - ${it.name}" }.join('\n')}

            Please ensure only one set of index files matches the reference basename pattern.
            """.stripIndent()
        }

        // check if all expected files are present
        def index_basename = basenames[0]

        // println "DEBUG: Using index basename: ${index_basename}"

        def expected_files = required_extensions.collect { ext ->
            file("${index_dir}/${index_basename}.${ext}")
        }
        def missing_files = expected_files.findAll { !it.exists() }
        def existing_files = expected_files.findAll { it.exists() }

        // println "DEBUG: Expected files:"
        //     expected_files.each { println "DEBUG:   - ${it.name} (exists: ${it.exists()})" }

        if (!missing_files.isEmpty()) {
            error """
            Incomplete BWA index for '${index_basename}' in ${index_dir}

            Missing files:
            ${missing_files.collect { "  - ${it.name}" }.join('\n')}

            Found files:
            ${existing_files.collect { "  - ${it.name}" }.join('\n')}

            A complete BWA index requires all 5 files: ${index_basename}.{amb,ann,bwt,pac,sa}
            Please rebuild the index using: bwa index ${ref}
            """.stripIndent()
        }

        // Create the channel with validated index directory (value channel allows singleton to be reused)
        ch_bwa_index = channel.value([
            [ id: index_basename ],
            index_dir
        ])

        // Alternative option using channel.of, requires collect() to create value channel, but the type changes to a list => not preferred
        // ch_bwa_index = channel.of(
        //     tuple(
        //         [ id: file(params.reference_fasta).simpleName ],
        //         file(params.reference_index, checkIfExists: true)
        //     )
        // ).collect()
        // ch_bwa_index = channel.fromPath(params.reference_fasta)
        //     .map( { ref ->  tuple([ id: ref.simpleName ], file(params.reference_index, checkIfExists: true)) } ).collect()
    }

    //
    // bwa alignment to reference genome
    //
    if (!skip_alignment) {
        sort_bam = true
        BWA_MEM (
            ch_fastq,
            ch_bwa_index,
            ch_ref_fasta,
            sort_bam        // markduplicates expects coordinate (or query) sorted input
        )
        ch_bam = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())

    // TODO: optionally enable CRAM output

    // TODO: when to combine runs/lanes from the same sample/library? how can stalling be avoided?
    // old sarek approach: https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L272
    // new sarek approach: https://github.com/nf-core/sarek/blob/20f41d1ce8b7ba296ee22adc71fe2da2ebcae93f/subworkflows/local/fastq_preprocess_gatk/main.nf#L163

        //
         // Mark duplicates using Picard and merge reads derived from same sample using RG values
        //

        // Group multi-lane/library samples by creating tuple with sample as grouping key
        ch_bam_grouped = ch_bam
            // remove unique distinguishing fields in meta data per sample
            .map { meta, bam ->
                [ meta - meta.subMap('id', 'read_group', 'single_end') + [ data_type: 'bam' ], bam ]
            }
            .groupTuple()
            .map { meta, bam -> [ [ id: meta.sample ] +  meta, bam ] }  // add new (sample-level) id

        GATK4_MARKDUPLICATES(
            ch_bam_grouped,
            ch_ref_fasta.map{ _meta, fasta -> [ fasta ] }.first(),
            ch_fai.map{ _meta, fai -> fai }.first()
        )
        ch_multiqc_files = ch_multiqc_files.mix(GATK4_MARKDUPLICATES.out.metrics.collect{it[1]})
        ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)
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
    //     ch_fai
    // )

    // TODO     SAMTOOLS_FAIDX.out.fai.map{ _meta, fai -> fai }.first() vs .collect()
    // To obtain a channel with just the file, which is not consumed by a task, a value channel is needed.
    // both first and collect create value channels, but first is clearer and creates a one item singleton channel, whereas collect creates a list.

    //
    // samtools sort and index bam files and collect alignment stats
    //

    // TODO: move into subworkflow?

        SAMTOOLS_SORT_MARKDUP(ch_bam_markdup, ch_ref_fasta, 'bai')
        ch_bam_markdup_sort = SAMTOOLS_SORT_MARKDUP.out.bam

        SAMTOOLS_INDEX(ch_bam_markdup_sort)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        ch_bam_bai = ch_bam_markdup_sort.join(SAMTOOLS_INDEX.out.bai)

        SAMTOOLS_STATS(ch_bam_bai, ch_ref_fasta)
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect{it[1]})

        SAMTOOLS_FLAGSTAT(ch_bam_bai)
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_FLAGSTAT.out.flagstat.collect{it[1]})
        ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)

        SAMTOOLS_IDXSTATS(ch_bam_bai)
        ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_IDXSTATS.out.idxstats.collect{it[1]})
        ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

        //
        // Module: mosdepth coverage statistics
        //

        ch_bam_bai_bed = ch_bam_bai.combine(ch_ref_bed.map { _meta, bed -> bed })
        MOSDEPTH(ch_bam_bai_bed, ch_ref_fasta)
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.per_base_bed.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_bed.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.quantized_bed.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.thresholds_bed.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.global_txt.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.regions_txt.collect{it[1]})
        ch_multiqc_files = ch_multiqc_files.mix(MOSDEPTH.out.summary_txt.collect{it[1]})
    }

    //
    // Call variants using GATK
    //

    if (!skip_variantcalling && !skip_alignment) {

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
        // VCF files (M interval VCFs, each with all N samples)
        //     ↓ PICARD_MERGEVCFS
        // Final VCF (1 file with all N samples, all M intervals)

        // Create GATK dictionary
        GATK4_CREATESEQUENCEDICTIONARY(ch_ref_fasta)
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
        ch_ref_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

        // Convert BED to GATK IntervalList
        GATK4_BEDTOINTERVALLIST(ch_ref_bed, ch_ref_dict)
        ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

        // TODO: add simpler strategy that just splits bed file into 1 task per contig + add more meta data about contig name

        // Split chrom/contigs into separate interval_list files for scatter-gather parallel processing
        GATK4_INTERVALLISTTOOLS(GATK4_BEDTOINTERVALLIST.out.interval_list)
        ch_versions = ch_versions.mix(GATK4_INTERVALLISTTOOLS.out.versions)
        ch_intervals = GATK4_INTERVALLISTTOOLS.out.interval_list
            // [ [interval_genome_meta.id], [interval_list, interval_list, ...] ] (single element)
            .transpose()
            // [ [interval_genome_meta.id], interval_list ] (multiple elements)

        ch_bam_bai_intervals = ch_bam_bai
            .combine(ch_intervals)
            .map { bam_meta, bam, bai, interval_genome_meta, interval ->
                def combined_meta = bam_meta + [
                    genome_id: interval_genome_meta.id,
                    interval_name: interval.simpleName
                ]
                [ combined_meta, bam, bai, interval ]
            }
            // [ [combined_meta.id, combined_meta.sample, combined_meta.num_entries, combined_meta.multiple_lanes, combined_meta.multiple_libraries, combined_meta.data_type, combined_meta.genome_id, combined_meta.interval_name], bam, bai, interval_list, [] ]

        // TODO: check if bed or interval_list is preferred

        // Alternative approach to generate scattered interval
        // BIN_INTERVALS(ch_ref_bed, 1000000, [])
        // TODO sarek alternative approach: https://github.com/nf-core/sarek/blob/master/modules/local/create_intervals_bed/main.nf
        // https://github.com/nf-core/sarek/blob/master/subworkflows/local/prepare_intervals/main.nf
        // https://github.com/nf-core/sarek/blob/master/subworkflows/local/bam_variant_calling_haplotypecaller/main.nf

        // Process each sample x interval combination in parallel with HaplotpeCaller
        // Group GVCFs for all samples per interval for GenomicsDBImport
        // Process each interval in parallel with GenotypeGVCFs
        // Merge intervals

        //
        // Module: GATK4_HAPLOTYPECALLER (sample x interval)
        // Call variants per sample per interval in parallel
        //
        GATK4_HAPLOTYPECALLER(
            ch_bam_bai_intervals.map{ combined_meta, bam, bai, interval ->
                [ combined_meta, bam, bai, interval, [] ] },
            ch_ref_fasta,
            ch_fai,
            ch_ref_dict,
            [[:], []],  // dbsnp (optional)
            [[:], []]   // dbsnp_tbi (optional)
        )
        ch_versions = ch_versions.mix(GATK4_HAPLOTYPECALLER.out.versions)
        ch_gvcf = GATK4_HAPLOTYPECALLER.out.vcf      // [[meta], vcf.gz]
        ch_gvcf_tbi = GATK4_HAPLOTYPECALLER.out.tbi  // [[meta], tbi]

        // Collect all samples for each interval for genomicsDBimport
        ch_gvcf_by_interval = ch_gvcf
            // [ [meta.id, meta.sample, meta.num_entries, meta.multiple_lanes, meta.multiple_libraries, meta.data_type, meta.genome_id, meta.interval_name], gvcf ]
            .join(ch_gvcf_tbi, by: 0)                   // Join gvcf and tbi by meta
            // [ meta, gvcf, tbi ]
            .map { meta, gvcf, tbi ->
                def interval_key = "${meta.genome_id}_${meta.interval_name}"
                [ interval_key, meta, gvcf, tbi ]
            }
            // [interval_key, meta, gvcf, tbi]
            .groupTuple(by: 0)                          // Group by interval
            // [ interval_key, [meta, meta, ...], [gvcf, gvcf, ...], [tbi, tbi, ...] ]
            .map { interval_key, metas, gvcfs, tbis ->
                def interval_meta = [                   // Create new meta without sample info
                    id: interval_key,
                    genome_id: metas[0].genome_id,
                    interval_name: metas[0].interval_name
                ]
                [ interval_meta, gvcfs, tbis ]
            }
            // [ interval_meta, [gvcfs], [tbis] ]

        // Rejoin with the actual interval file
        ch_genomicsdb_input = ch_gvcf_by_interval   // [ interval_meta, [gvcfs], [tbis] ]           (1 per sample)
            .combine(ch_intervals)                  // [ [interval_genome_meta.id], interval_list ] (1 per interval)
            // [ interval_meta, [gvcfs], [tbis], [interval_genome_meta.id], interval_list]          (sample x interval)
            .map { interval_meta, gvcfs, tbis, _interval_meta, interval_file ->
                // Match by interval name
                if (interval_file.simpleName == interval_meta.interval_name) {
                    return [ interval_meta, gvcfs, tbis, interval_file, [], []]
                }
            }
            .filter { it != null }  // Remove non-matching combinations
            // TODO: is this needed?

        //
        // Module: GATK4_GENOMICSDBIMPORT (intervals)
        // Consolidate GVCFs per interval across all samples
        //
        GATK4_GENOMICSDBIMPORT(
            ch_genomicsdb_input,
            false,  // run_intlist
            false,  // run_updatewspace
            false   // input_map
        )
        ch_versions = ch_versions.mix(GATK4_GENOMICSDBIMPORT.out.versions)
        ch_genomicsdb = GATK4_GENOMICSDBIMPORT.out.genomicsdb  // [[meta], genomicsdb_dir]

        // Prepare channels for GenotypeGVCFs by matching genomicsdb with corresponding interval file
        ch_genotype_input = ch_genomicsdb   // (1 per interval)
            .combine(ch_intervals)          // (1 per interval)
            // [ meta, genomicsdb_dir, interval_meta, interval_list ] (interval x genomicsb_interval_dirs => contains mismatched elements
            .map { meta, genomicsdb, _interval_meta, interval_file ->
                // Match by interval name
                if (interval_file.simpleName == meta.interval_name) {
                    return [ meta, genomicsdb, [], interval_file, [] ]
                }
            }
            .filter { it != null }

        //
        // Module: GATK4_GENOTYPEGVCFS (intervals)
        // Perform joint genotyping per interval
        //
        GATK4_GENOTYPEGVCFS(
            ch_genotype_input,
            ch_ref_fasta,
            ch_fai,
            ch_ref_dict,
            [[:], []],  // dbsnp (optional)
            [[:], []]   // dbsnp_tbi (optional)
        )
        ch_versions = ch_versions.mix(GATK4_GENOTYPEGVCFS.out.versions)
        ch_vcf_by_interval = GATK4_GENOTYPEGVCFS.out.vcf  // [[meta], vcf.gz]
        ch_vcf_tbi_by_interval = GATK4_GENOTYPEGVCFS.out.tbi

        // Prepare for VCF merging: collect all interval VCFs
        // Group by genome_id to merge all intervals for that genome
        ch_vcfs_to_merge = ch_vcf_by_interval
            // [ meta, vcf ]    (1 per interval)
            .join(ch_vcf_tbi_by_interval, by: 0)
            // [ meta, vcf, tbi ]    (1 per interval)
            .map { meta, vcf, tbi ->
                [ meta.genome_id, vcf, tbi ]
            }
            .groupTuple(by: 0)
            // [ meta.genome_id, [vcfs],  [tbis] ] (1 element with all interval files)
            .map { genome_id, vcfs, _tbis ->
                def final_meta = [ id: genome_id ]
                [ final_meta, vcfs ]
            }

        //
        // Module: GATK4_MERGEVCFS
        // Merge interval VCFs into single cohort VCF
        // Order is handled automatically by the sequence dictionary
        // (unlike bcftools, see https://nf-co.re/subworkflows/vcf_gather_bcftools)
        //
        GATK4_MERGEVCFS(
            ch_vcfs_to_merge,
            ch_ref_dict
        )
        ch_versions = ch_versions.mix(GATK4_MERGEVCFS.out.versions)
        ch_final_vcf = GATK4_MERGEVCFS.out.vcf      // [[meta], merged.vcf.gz]
        ch_final_vcf_tbi = GATK4_MERGEVCFS.out.tbi  // [[meta], merged.vcf.gz.tbi]

        //
        // Module: TABIX_TABIX
        // Index the final merged VCF
        //
        // Or use bcftools? https://nf-co.re/modules/bcftools_index/
        TABIX_TABIX(ch_final_vcf)
        ch_final_vcf_tbi = TABIX_TABIX.out.index
        // TODO not required?
    }


    // TODO add BQSR

    // TODO: add VQSR + hard filtering

    // TODO: fastq-screen -> see separate branch

    // TODO: filter vcf

    //
    // snpEff annotation
    //

    if (!skip_annotation && !skip_variantcalling && !skip_alignment) {
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
    }

    // TODO:     gatk VariantsToTable -V "${ann_dir}/${species}/combined.filter_added.ann.vcf" -F CHROM -F POS -F TYPE -GF GT -O "${ann_dir}/${species}/combined.filter_added.table"

    // TODO: bcf normalizes

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

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'plasmovar_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        channel.fromPath(params.multiqc_config, checkIfExists: true) :
        channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
