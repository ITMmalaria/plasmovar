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
include { BWA_INDEX                              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                                } from '../modules/nf-core/bwa/mem/main'
include { paramsSummaryMap                       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_plasmovar_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DEFINE ADDITIONAL FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


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

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // TODO create channel holding reference genome
    // Note: bwa_index expects meta channel, whereas bbsplit_indexer just needs a path
    // ch_reference = Channel.value(file(params.reference))
    // ch_reference = Channel.fromPath(params.reference).map{ref -> tuple (ref.simpleName, ref)}
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

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_samplesheet
    )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})  // MultiQC's fastqc module requires the zip output - https://multiqc.info/modules/fastqc/
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

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
        ch_adapter_fasta = params.fastp_adapter_fasta ? Channel.fromPath(param, checkIfExists: true).collect() : []
        FASTP (
            ch_samplesheet, // channel: [ val(meta), [ reads ] ]
            ch_adapter_fasta,
            params.fastp_save_trimmed,
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
    // Host read filtering / host sequence contamination removal / host decontamination
    //

    // TODO: add BWA decontamination option
    // TODO: move to subworkflow
    // TODO: use file with list of fasta paths and names? samplesheet could contain species
    // cf. samplesheet could mention which species to map against in bbsplit (human) and fastq screen
    // TODO: examples:
    // https://github.com/nf-core/eager/blob/dev/modules/local/host_removal.nf
    // https://github.com/nf-core/taxprofiler/blob/1.1.7/subworkflows/local/shortread_hostremoval.nf
    if (!params.skip_hostremoval) {
        // use provided BBSplit index if supplied or generate from scratch otherwise
        if (!params.hostremoval_bbsplit_index) {
            // Prepare channel with list of reference genomes to filter reads against.
            // Expected format for BBMAP_BBSPLIT module is:
            //      tuple val(other_ref_names), path (other_ref_paths)
            //      [['name'], [/path/to/fast.gz]]
            Channel.from( [
                [params.hostremoval_bbsplit_reference_name,
                params.hostremoval_reference]
            ] )
                .collect{ id, fasta -> [ [id], [file(fasta, checkIfExists: true)] ] }
                .set { ch_bbsplit_other_refs }


            // create bbsplit index for filtering
            BBMAP_BBSPLIT_INDEXER (
                [ [:], [] ],
                [],
                Channel.value(file(params.reference, checkIfExists: true)),
                ch_bbsplit_other_refs,
                true
            )
            ch_bbsplit_index = BBMAP_BBSPLIT_INDEXER.out.index
            ch_versions = ch_versions.mix(BBMAP_BBSPLIT_INDEXER.out.versions)
            // bbsplit.sh -Xmx6000M ref_primary="/path/to/primary_genome.fasta"  ref_human="/path/to/contaminant_genome.fa.gz" path=bbsplit_index_output threads=4
        } else {
            // Index needs to be the directory `genome/index/bbsplit` which contains a ref subdir,
            // which in turn contains an index and genome subdir.
            Sytem.println("Using pre-supplied reference fasta for host removal")
            ch_bbsplit_index = Channel.value(file(params.hostremoval_bbsplit_index, checkIfExists: true))
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
    } else {
        ch_reads_for_alignment = ch_reads_for_hostremoval
    }

    //
    // MODULE: Run bwa index
    // Index reference genome
    //
    if (!params.reference_index) {
        BWA_INDEX (
            Channel.fromPath(params.reference).map{ref -> tuple (ref.simpleName, ref)}
        )
        ch_bwa_index = BWA_INDEX.out.index.collect()
        // TODO: use .versions.first() or just .versions?
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions.first())
    } else {
        ch_bwa_index = Channel.value(file(params.reference_index, checkIfExists: true))
    }
    // TODO conditional exit or only use skips?
    if (params.only_build_reference) {
        System.println("conditional check only build reference")
        System.exit(0)
    }

    //
    // MODULE: Run bwa mem
    // Alignment to reference genome
    //
    sort_bam = true
    BWA_MEM (
        ch_reads_for_alignment,
        ch_bwa_index,
        // ch_bwa_index.map{ it -> [ [ id:'index' ], it ] },
        [[id:'no_fasta'], []],
        // [[],[]],
        // sort
        sort_bam
    )
    ch_versions = ch_versions.mix(BWA_MEM.out.versions.first())
    // logs?

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
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

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
