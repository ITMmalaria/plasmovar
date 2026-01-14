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

include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { FASTP                                  } from '../modules/nf-core/fastp/main'
include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_INDEXER } from '../modules/nf-core/bbmap/bbsplit/main'
include { BBMAP_BBSPLIT as BBMAP_BBSPLIT_MAPPER  } from '../modules/nf-core/bbmap/bbsplit/main'
include { DEACON_INDEX                           } from '../modules/nf-core/deacon/index/main'
include { DEACON_INDEX_DIFF                      } from '../modules/local/deacon/diff/main'
include { DEACON_FILTER                          } from '../modules/nf-core/deacon/filter/main'
include { BWA_INDEX                              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                                } from '../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_FAIDX                         } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_INDEX                         } from '../modules/nf-core/samtools/index/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT_MARKDUP } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_STATS                         } from '../modules/nf-core/samtools/stats/main'
include { SAMTOOLS_FLAGSTAT                      } from '../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_IDXSTATS                      } from '../modules/nf-core/samtools/idxstats/main'
include { GATK4_MARKDUPLICATES                   } from '../modules/nf-core/gatk4/markduplicates/main'
include { CREATE_INTERVALS_BED                   } from '../modules/local/create_intervals_bed/main'
include { MOSDEPTH                               } from '../modules/nf-core/mosdepth/main'
include { GATK4_CREATESEQUENCEDICTIONARY         } from '../modules/nf-core/gatk4/createsequencedictionary/main'
include { GATK4_BEDTOINTERVALLIST                } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_INTERVALLISTTOOLS                } from '../modules/nf-core/gatk4/intervallisttools/main'
include { BIN_INTERVALS                          } from '../modules/local/bin_intervals/main'
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
// TODO: add unit test
def addReadgroupToMeta(meta, files) {
    // Note: needs to be run on initial file path before processing,
    // otherwise the path in the work dir will not be found
    def flowcell = flowcellLaneFromFastq(files[0])
    if ( !meta.single_end && flowcell != flowcellLaneFromFastq(files[1]) ){
        error("Flowcell ID does not match for paired reads of sample ${meta.id} - ${files}")
    }

    // Define RG values, only adding non-empty values
    // This avoids problems with Picard and the need for the `VALIDATION_STRINGENCY LENIENT` option,
    // since there will not be any empty fields like:
    // @RG	ID:106264-002-079.22NY35LT3.L007	CN:	PU:22NY35LT3.L007	SM:106264-002-079	LB:106264-002-079	PL:ILLUMINA

    // If we cannot read the flowcell ID from the fastq file, then we don't use it
    def RG_ID = flowcell ? "${meta.sample}.${flowcell}.${meta.lane}" : "${meta.sample}.${meta.lane}"

    // For LB, check if library is defined in input samplesheet (trim handles case of empty string)
    def RG_LB = meta.library?.trim() ? "${meta.sample}.${meta.library}" : "${meta.sample}"

    def RG_map = [
        ID: RG_ID,
        SM: meta.sample,
        LB: RG_LB,
        PL: params.seq_platform,                            // defaults to ILLUMINA
        CN: params.seq_center?.trim() ?: null,              // most likely blank
        PU: flowcell ? "${flowcell}.${meta.lane}" : null,   // used by BQSR if present, otherwise it defaults to using ID
    ]

    // Build RG string: filter out empty fields and format as RG string
    def RG_fields = RG_map
            // example input: [ID: "sample1", SM: "sample1", CN: null, PU: "flowcell.L001"]
        .findAll { _key, value ->
            value != null && value != ''    // Explicit: keep only non-empty values
        }   // returns map: [ID: "sample1", SM: "sample1", PU: "flowcell.L001"]
        .collect { key, value ->
            "${key}:${value}"               // Format each entry as "KEY:VALUE"
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

    // TODO: set correct skip options depending on which "*_only" options were enabled
    // params.only_build_reference -> skip_trimming, skip_hostremoval, skip_alignment, skip_variantcalling

    // TODO create channel holding reference genome
    // Note: bwa_index expects meta channel, whereas bbsplit_indexer just needs a path
    // ch_reference = channel.value(file(params.reference))
    // ch_reference = channel.fromPath(params.reference).map{ref -> tuple (ref.simpleName, ref)}
    // reference = [[ id:'reference', primary:true ], file(params.reference)]



    // TODO: allow either a single fasta reference file to be supplied or multiple ones
    // via samplesheet column or via extra file listing species name + reference path
    // optionally providing pre-built indices (bwa, samtools, gatk, etc.)
    // for reference sheet example: https://nf-co.re/eager/dev/docs/usage#reference-input
    // also requires changes to bbsplit (needs both host and main reference)

    // TODO: integrate all modules/subworkflows into multiQC: https://nf-co.re/docs/contributing/tutorials/adding_modules_to_pipelines

    // TODO: bundle fastq preprocessing steps into subworkflow

    // TODO: https://nf-co.re/modules/fastqscreen_fastqscreen

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    // TODO: see nf-core/rnaseq and decide when/where concatenation should happen

    // TODO: see description of from samplesheet input validation plugin:
    //https://github.com/nextflow-io/nf-validation/blob/750a56d02ce902508eb7777188b034d0b8f3435c/docs/samplesheets/examples.md
    // TODO: also check structure of meta map
    // TODO: describe columns as done here https://nf-co.re/sarek/3.4.4/docs/usage/


    // Add read groups to meta
    if (!params.only_build_reference) {
        ch_samplesheet = ch_samplesheet.map { meta, fastqs -> addReadgroupToMeta(meta, fastqs) }
    }

    if (!params.skip_qc) {

        //
        // MODULE: Run FastQC
        //
        FASTQC (
            ch_samplesheet
        )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})  // MultiQC's fastqc module requires the zip output - https://multiqc.info/modules/fastqc/
    }

    //
    // MODULE: FASTP - trim reads using fastp
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
    if (!params.skip_trimming) {
        // create expected input for fastp module using read adapter file path from input parameters
        // channel: [ val(meta), [ path(reads) ], path(adapters) ]
        ch_samplesheet
            .map { meta, reads -> [ meta, reads, params.fastp_adapter_fasta ?: [] ] }
            .set { ch_fastp_input }

        FASTP (
            ch_fastp_input,
            false,  // discard_trimmed_pass - Specify true to not write any reads that pass trimming thresholds. This can be used to use fastp for the output report only. Previously set to `!params.fastp_save_trimmed`
            params.fastp_save_trimmed_fail,
            params.fastp_save_merged
        )
        ch_reads_for_hostremoval = FASTP.out.reads
        ch_versions = ch_versions.mix(FASTP.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTP.out.json.collect{it[1]})  // MultiQC's fastp module relies only on the json output - https://multiqc.info/modules/fastp/
    // TODO: re-run fastQC after fastp, or just rely on subworkflow mentioned above
    } else {
        ch_reads_for_hostremoval = ch_samplesheet
    }
    // TODO: use more descriptive name in else clause to avoid wrongly named channels
    // sarek uses "reads_for_nexttool"

    //
    // Host removal / depletion / read filtering / host sequence contamination removal / host decontamination
    //

    // TODO: add BWA decontamination option
    // TODO: move to subworkflow?
    // TODO: use file with list of fasta reference paths and names?
    // samplesheet could contain species
    // cf. samplesheet could mention which species to map against in bbsplit (human) and fastq screen
    // TODO: examples:
    // https://github.com/nf-core/eager/blob/dev/modules/local/host_removal.nf
    // https://github.com/nf-core/taxprofiler/blob/1.1.7/subworkflows/local/shortread_hostremoval.nf

    if (!params.skip_hostremoval) {

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
                    channel.value(file(params.reference, checkIfExists: true)),
                    ch_bbsplit_other_refs,
                    true
                )
                ch_bbsplit_index = BBMAP_BBSPLIT_INDEXER.out.index
                ch_versions = ch_versions.mix(BBMAP_BBSPLIT_INDEXER.out.versions.first())
                // bbsplit.sh -Xmx6000M ref_primary="/path/to/primary_genome.fasta"  ref_human="/path/to/contaminant_genome.fa.gz" path=bbsplit_index_output threads=4
            } else {
                // Index needs to be the directory `genome/index/bbsplit` which contains a ref subdir,
                // which in turn contains an index and genome subdir.
                System.println("Using pre-supplied reference fasta for host removal")
                ch_bbsplit_index = channel.value(file(params.hostremoval_bbsplit_index, checkIfExists: true))
            }

            // run bbsplit in map mode
            BBMAP_BBSPLIT_MAPPER (
                ch_reads_for_hostremoval,
                ch_bbsplit_index,
                [],
                [ [], [] ],
                false
            )
            ch_reads_for_alignment = BBMAP_BBSPLIT_MAPPER.out.primary_fastq
            ch_versions = ch_versions.mix(BBMAP_BBSPLIT_MAPPER.out.versions.first())
            // ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBSPLIT_MAPPER.out.stats.collect{it[1]})
            // TODO multiqc bbsplit not showing up due to bug https://github.com/MultiQC/MultiQC/pull/1513
        }
        // use Deacon for host read removal
        // TODO: add prebuilt pangenome index https://github.com/bede/deacon?tab=readme-ov-file#prebuilt-indexes
        else if (params.hostremoval_method == "deacon") {

            // create Deacon index if not provided
            if (!params.hostremoval_deacon_index) {

                def deacon_fasta = file(params.hostremoval_reference)
                ch_deacon_fasta = channel.of(
                    tuple(
                        [ id: deacon_fasta.simpleName ],
                        file(deacon_fasta)
                    )
                )

                DEACON_INDEX (
                    ch_deacon_fasta
                )

                ch_reads_for_hostremoval
                    .combine(DEACON_INDEX.out.index)
                    .map { meta, reads, _meta_index, index -> [ meta, index, reads ] }
                    .set { ch_deacon_input }

                ch_versions = ch_versions.mix(DEACON_INDEX.out.versions.first())

            } else {
                // retrieve index from input parameters
                ch_reads_for_hostremoval
                    .map { meta, reads -> [ meta, params.hostremoval_deacon_index ?: [] ,  reads] }
                    .set { ch_deacon_input }
            }

            if ( params.hostremoval_deacon_diff_fasta ) {
                def deacon_diff_fasta = file(params.hostremoval_deacon_diff_fasta)

                // returns too many elements
                // ch_deacon_diff = ch_deacon_input
                //     .map { meta, index, reads ->
                //         tuple(
                //             [ id: deacon_diff_fasta.baseName ],
                //             index,
                //             deacon_diff_fasta
                //         )
                //     }

                ch_deacon_index = ch_deacon_input
                    // .map { tuple -> tuple[1] }  // Extract the index from each tuple (index is the second element)
                    .map { _meta, index, _reads -> [index] }
                    .unique()

                ch_deacon_diff = ch_deacon_index
                    .map { index ->
                        tuple(
                            [ id: "${index.simpleName}_diff_${deacon_diff_fasta.simpleName}" ],
                            index,
                            deacon_diff_fasta
                        )
                    }

                DEACON_INDEX_DIFF(ch_deacon_diff)

                ch_deacon_input = ch_deacon_input
                    .combine(DEACON_INDEX_DIFF.out.index)
                    .map { meta, _index, reads, _meta_diff, diff_index -> [ meta, diff_index, reads] }
                ch_versions = ch_versions.mix(DEACON_INDEX_DIFF.out.versions.first())
            }

            // filter reads using deacon against host index
            DEACON_FILTER(ch_deacon_input)
            ch_reads_for_alignment = DEACON_FILTER.out.fastq_filtered
            ch_versions = ch_versions.mix(DEACON_FILTER.out.versions.first())
        }
    } else {
        // skip host removal and continue unmodified read channel for alignment
        ch_reads_for_alignment = ch_reads_for_hostremoval
    }

    //
    // MODULE: Run bwa index
    // Index reference genome if index is not already provided
    //

    // Create channel containing the reference fasta
    // ch_ref_fasta = channel.fromPath(params.reference)
    //     .map { ref -> tuple([ id: ref.simpleName ], ref) }
    //     .collect()   // or .first()
    def ref = file(params.reference, checkIfExists: true)
    def ref_basename = ref.simpleName
    ch_ref_fasta = channel.value([
        [ id: ref_basename ],
        ref
    ])

    if (!params.reference_index) {
        // Construct bwa index for the reference fasta if it is not supplied by the user
        BWA_INDEX (ch_ref_fasta)
        ch_bwa_index = BWA_INDEX.out.index.collect()    // collect() is required to create a re-usable value channel, otherwise it will only contain a single element which won't be emitted for each of the sample reads in ch_reads_for_alignment
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

        // Create the channel with validated index directory
        ch_bwa_index = channel.value([
            [ id: index_basename ],
            index_dir
        ])

        // Alternative option using channel.of, requires collect() to create value channel
        // ch_bwa_index = channel.of(
        //     tuple(
        //         [ id: file(params.reference).simpleName ],
        //         file(params.reference_index, checkIfExists: true)
        //     )
        // ).collect()
        // ch_bwa_index = channel.fromPath(params.reference)
        //     .map( { ref ->  tuple([ id: ref.simpleName ], file(params.reference_index, checkIfExists: true)) } )
    }

    //
    // MODULE: Run bwa mem
    // Alignment to reference genome
    //
    if (!params.skip_alignment) {
        sort_bam = true
        BWA_MEM (
            ch_reads_for_alignment,
            ch_bwa_index,
            ch_ref_fasta,
            sort_bam        // markduplicates expects coordinate (or query) sorted input
        )
        ch_bam = BWA_MEM.out.bam
        ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    }
    // TODO: optionally enable CRAM output

    // TODO: when to combine runs/lanes from the same sample/library? how can stalling be avoided?
    // old sarek approach: https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L272
    // new sarek approach: https://github.com/nf-core/sarek/blob/20f41d1ce8b7ba296ee22adc71fe2da2ebcae93f/subworkflows/local/fastq_preprocess_gatk/main.nf#L163

    //
    // MODULE: Run picard markduplicates
    // Mark duplicates using Picard and merge reads derived from same sample using RG values
    //

    // Index reference fasta as .fai
    SAMTOOLS_FAIDX(
        ch_ref_fasta,
        [[:], []],
        []
    )
    ch_versions = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)
    ch_fai      = SAMTOOLS_FAIDX.out.fai    // [[id:Pf3D7_01_v3], /path/to/ref.fa.gz.fai]

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
    ch_versions = ch_versions.mix(GATK4_MARKDUPLICATES.out.versions)
    ch_bam_markdup = GATK4_MARKDUPLICATES.out.bam

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
    // Module: samtools sorting, indexing and stats collection
    //

    // TODO: move into subworkflow?

    SAMTOOLS_SORT_MARKDUP(ch_bam_markdup, ch_ref_fasta, 'bai')
    ch_bam_markdup_sort = SAMTOOLS_SORT_MARKDUP.out.bam

    SAMTOOLS_INDEX(ch_bam_markdup_sort)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
    ch_bam_bai = ch_bam_markdup_sort.join(SAMTOOLS_INDEX.out.bai)

    SAMTOOLS_STATS(ch_bam_bai, ch_ref_fasta)
    SAMTOOLS_FLAGSTAT(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions)
    SAMTOOLS_IDXSTATS(ch_bam_bai)
    ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions)

    //
    // Module: mosdepth coverage statistics
    //

    if (params.intervals_bed) {
        bed_file = file(params.reference, checkIfExists: true)
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
    }

    ch_bam_bai_bed = ch_bam_bai.combine(ch_ref_bed.map { _meta, bed -> bed })
    MOSDEPTH(ch_bam_bai_bed, ch_ref_fasta)

    //
    // Module: GATK variant calling
    //

    // Create GATK dictionary
    GATK4_CREATESEQUENCEDICTIONARY(ch_ref_fasta)
    ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    ch_ref_dict = GATK4_CREATESEQUENCEDICTIONARY.out.dict

    // Convert BED to GATK IntervalList
    GATK4_BEDTOINTERVALLIST(ch_ref_bed, ch_ref_dict)
    ch_versions = ch_versions.mix(GATK4_BEDTOINTERVALLIST.out.versions)

    // Split chrom/contigs into separate interval_list files for scatter-gather parallel processing
    GATK4_INTERVALLISTTOOLS(GATK4_BEDTOINTERVALLIST.out.interval_list)
    ch_versions = ch_versions.mix(GATK4_INTERVALLISTTOOLS.out.versions)
    ch_intervals = GATK4_INTERVALLISTTOOLS.out.interval_list

    // Alternative approach to generate scattered interval
    // BIN_INTERVALS(ch_ref_bed, 1000000, [])
    // TODO sarek alternative approach: https://github.com/nf-core/sarek/blob/master/modules/local/create_intervals_bed/main.nf
    // https://github.com/nf-core/sarek/blob/master/subworkflows/local/prepare_intervals/main.nf
    // https://github.com/nf-core/sarek/blob/master/subworkflows/local/bam_variant_calling_haplotypecaller/main.nf

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
