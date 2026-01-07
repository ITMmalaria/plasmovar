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
        // [[id:sample, lane:1, library:null], /path/to/sample_L001_R1.fastq.gz, /path/to/sample_L001_R2.fastq.gz]
        // [[sample:sample_1, lane:1], /path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]
        // extract sample id (can be shared by multiple fastq (pairs)), so we can group by it
        .map { meta, fastq_1, fastq_2 ->
            // auto-detect lane if not provided
            if (params.auto_detect_lanes && !meta.lane && meta.id) {
                def detected_lane = extractLaneFromFilename(fastq_1.toString())
                if (detected_lane) {
                    meta = meta + [lane: detected_lane]
                    log.warn("Auto-detected lane ${detected_lane} for sample ${meta.id} with fastq file(s) ${fastq_1} ${fastq_2}")
                }
            }
            // handle single and paired-end fastq files
            if (!fastq_2) {
                return [ meta.id, [ meta + [ single_end: true, sample: "${meta.id}" ], [ fastq_1 ] ] ]
            } else {
                return [ meta.id, [ meta + [ single_end: false, sample: "${meta.id}" ], [ fastq_1, fastq_2 ] ] ]
            }
            // TODO compare with nf-core template approach, which uses fewer nested lists + flatten.
        }
        .tap{ ch_samples } // save sample-wise channel for later re-use
        // [sample_1, [[id:sample_1, lane:1, library: null, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]]
        // group by sample id
        .groupTuple()
        // [sample_1, [[[id:sample_1, lane:1, library:null, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]], [[id:sample_1, lane:2, library:null, single_end:false, sample:sample_1], [/path/to/sample_1_L002_R1_001.fastq.gz, /path/to/sample_1_L002_R2_001.fastq.gz]]]]
        .map { validateInputSamplesheet(it) }
        // calculate number of fastq (pairs) per sample
        .map { sample, ch_items ->
            def metas = ch_items.collect { it[0] }
            def num_entries = ch_items.size()
            def has_multiple_lanes = metas.collect { it.lane }.findAll { it != null }.unique().size() > 1
            def has_multiple_libraries = metas.collect { it.library }.findAll { it != null }.unique().size() > 1
            return [ sample, num_entries, has_multiple_lanes, has_multiple_libraries ]
        }
        // [sample_1, 2, true, false]
        .combine(ch_samples, by: 0)
        // [sample_1, 2, true, false, [[id:sample_1, lane:1, library:null, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]]
        // [sample_1, 2, true, false, [[id:sample_1, lane:2, library:null, single_end:false, sample:sample_1], [/path/to/sample_1_L002_R1_001.fastq.gz, /path/to/sample_1_L002_R2_001.fastq.gz]]]
        .map { _sample, num_entries, multi_lane, multi_lib, ch_item ->
            def ( meta, fastqs ) = ch_item

            // build unique identifier based on sample name, library and lane
            def id_elements = [meta.sample]
            if (meta.library) id_elements.add("LIB_${meta.library}")
            if (meta.lane) id_elements.add("LANE_${meta.lane}")

            meta = meta + [
                id: id_elements.join('-'),
                num_entries: num_entries.toInteger(),
                multiple_lanes: multi_lane,
                multiple_libraries: multi_lib
            ]

            return [ meta, fastqs ]
        }
        // [[id:sample_1-LANE_1, lane:1, library:null, single_end:false, sample:sample_1, num_entries:2, multiple_lanes:true, multiple_libraries:false], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]
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
    def (sample, items) = input
    def metas = items.collect { it[0] }

    // Only validate if there are multiple fastq (pair) entries for the same sample
    if (metas.size > 1) {
        // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
        def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
        if (!endedness_ok) {
            error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
        }

        def has_lane_info = metas.every { it.lane != null }
        def has_library_info = metas.every { it.library != null }

        // Check if samples have distinguishing info
        if (!has_lane_info && !has_library_info) {
            error """
                Sample '${sample}' has multiple entries but no lane or library information.
                Please provide either 'lane' or 'library' columns in your samplesheet to distinguish between files.
            """.stripIndent()
        }

        // Check if lane+library combinations are unique
        def identifiers = metas.collect { "${it.lane ?: 'NA'}_${it.library ?: 'NA'}" }
        if (identifiers.size() != identifiers.unique().size()) {
            error """
                Sample '${sample}' has duplicate lane/library combinations.
                Each sample must have a unique combination of lane and library identifiers.
            """.stripIndent()
        }

        // TODO: Can optionally be replaced by adding the following to schema_input.json:
        // "uniqueEntries": ["sample", "lane", "library"]
        // Should be placed just after the items declaration (i.e. level 0).
        // However, error message will be less detailed.
        // Also needs verification whether shared null/[] entries would trigger an error or not.

        // TODO: write unit tests for these checks
    }
    return input
}

//
// Attempt to extract lane info from fastq filename, as a fallback when this info is not provided
//
def extractLaneFromFilename(filename) {
    // [_\.] - ensures lane is preceded/followed by underscore or dot (not just anywhere)
    // (\d{1,3}) - capture group for the lane number between 1 and 3 digits long
    def patterns = [
        ~/.*[_\.]L(\d{1,3})[_\.].*/,            // _L001_, _L1_, .L01.
        ~/.*[_\.]lane[_\-]?(\d{1,3})[_\.].*/    // _lane1_ or _lane-01_
    ]

    return patterns.findResult { pattern ->
        def matcher = filename =~ pattern
        if (matcher) {
            def lane = matcher[0][1]    // matcher format = [sample_L001_R1.fastq.gz, 001]
            return lane.padLeft(3, '0') // pad lane to 3 digits for consistency
        }
        return null // return nothing if no pattern is found
    }

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
