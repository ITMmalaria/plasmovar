//
// Subworkflow with functionality specific to the pmoris/plasmovar pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsHelp                } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet
    help              // boolean: Display help message and exit
    help_full         // boolean: Show the full help message
    show_hidden       // boolean: Show hidden parameters in the help message

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        "",
        "",
        command
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Create channel from input file provided through params.input
    //
    // Adapted from: https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/subworkflows/local/samplesheet_to_channel/main.nf
    channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        // [[id:sample, lane:1], /path/to/sample_L001_R1.fastq.gz, /path/to/sample_L001_R2.fastq.gz]
        // [[sample:sample_1, lane:1], /path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]
        .map {
            meta, fastq_1, fastq_2 ->
                if (!fastq_2) {
                    return [ meta.id, [ meta + [ single_end:true, sample:"$meta.id" ], [ fastq_1 ] ] ]
                } else {
                    return [ meta.id, [ meta + [ single_end:false, sample:"$meta.id" ], [ fastq_1, fastq_2 ] ] ]
                }
        }.tap{ ch_samples } // save sample-wise channel
        // [sample_1, [[id:sample_1, lane:1, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]]
        .groupTuple()   // group by sample id -
        // [sample_1, [[[id:sample_1, lane:1, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]], [[id:sample_1, lane:2, single_end:false, sample:sample_1], [/path/to/sample_1_L002_R1_001.fastq.gz, /path/to/sample_1_L002_R2_001.fastq.gz]]]]
        // .view()
        .map { validateInputSamplesheet(it) }
        .map { sample, ch_items -> [ sample, ch_items.size() ] } // get number of lanes per sample
        .combine(ch_samples, by: 0) // for each entry add numLanes
        // [sample_1, 2, [[id:sample_1, lane:1, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]]
        .map {
            sample, num_lanes, ch_items ->
            ( meta, fastqs ) = ch_items
            if (meta.lane) {
                meta = meta + [id: "${meta.sample}-${meta.lane}".toString(), num_lanes: num_lanes.toInteger()]
            }
            // no need for else statement because id is already included from the start
            return [ meta, fastqs ]
        }
        // [[id:sample_1-1, lane:1, single_end:false, sample:sample_1, num_lanes:2], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]
        .set { ch_samplesheet }

        // TODO: check new nfcore/tools method:
        // .groupTuple()
        // .map { samplesheet ->
        //     validateInputSamplesheet(samplesheet)
        // }
        // .map {
        //     meta, fastqs ->
        //         return [ meta, fastqs.flatten() ]

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

// TODO: lane can be changed to "run" to be more general?

// TODO: check what happens if there is no lane info (or no laned runs): https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L287 (meta.size might be needed, but it is set to 1 by default in sarek for all entries unless something happens in align step?)

// TODO: check behaviour for single end data

// TODO: num_lanes can be removed if not used for checking output bams: https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/conf/modules/aligner.config#L51
// if so, entire block can be simplified by just creating meta.sample-meta.lane immediately

// TODO; try to create sample-lane id manual, then add num_lanes, check what size= is, add branching option for single end reads

// TODO: changed to sample, re-create id (combined sample-lane) because this is used by downstream modules

// TODO: remove first part because of redundancy, not necessary if there are no patients....although we do need the count of the number of lanes per sample
// TODO or do we?
// see readgroup construct: https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L949

// TODO: check branching options for multiple types of inputs https://github.com/nextflow-io/nf-validation/blob/750a56d02ce902508eb7777188b034d0b8f3435c/docs/samplesheets/examples.md

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    // Extract meta arrays - expected format:
    // [sample_1, [[[id:sample_1, lane:1, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]], [[id:sample_1, lane:2, single_end:false, sample:sample_1], [/path/to/sample_1_L002_R1_001.fastq.gz, /path/to/sample_1_L002_R2_001.fastq.gz]]]]
    def metas = input[1].collect{ it[0] }

    // Perform checks if there are multiple runs for the same sample
    if (metas.size > 1) {
        // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
        def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
        if (!endedness_ok) {
            error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
        }
        // and that lanes (if provided) are unique
        def lanes_ok = metas.collect{ it.lane }.unique().size == metas.size
        // def lanes_ok = metas.collect{ it.lane }.unique().size == fastqs.collect{ it[0] }.size   // equivalent to fastqs.size(), but not fastqs.size (=> latter gives length of individual fastq arrays in bag)
        if (!lanes_ok) {
            error("Please check input samplesheet -> Multiple runs of a sample must have a different lane or run number to differentiate between them: ${metas[0].id}")
        }
        // TODO: write unit tests for these checks
    }
    return input
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
