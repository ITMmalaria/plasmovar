//
// Subworkflow with functionality specific to the ITMmalaria/plasmovar pipeline
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

// TODO: check input vs params.input and monochrome_logs options. Might get fixed after nf-core tools update.

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

    def before_text = ""
    def after_text = ""
    if (monochrome_logs) {
        before_text = before_text.replaceAll(/\033\[[0-9;]*m/, '')
    }

    command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"

    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null,
        help,
        help_full,
        show_hidden,
        before_text,
        after_text,
        command,
        null
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
            // auto-detect lane if not provided - used for checking uniqueness and read group construction
            if (params.auto_detect_lanes && !meta.lane) {
                // TODO: add fallback option for extracting lane info from fastq identifiers

                def detected_lane = extractLaneFromFilename(fastq_1.toString())

                if (detected_lane) {
                    log.warn("Auto-detected lane ${detected_lane} for sample ${meta.id} with fastq file(s) ${fastq_1} ${fastq_2}")
                }
                else {
                    log.warn("Failed to extract lane info from FASTQ file name: ${fastq_1}")
                    if (! params.strict_lane_detection) {
                        log.warn("Proceeding without lane info for ${meta.id}, but this could lead to duplicate input warnings or incorrect results.")
                    } else {
                        log.error("Stopping pipeline since strict_lane_detection was enabled.")
                    }
                }

                // check if detected lanes match for paired reads
                if (fastq_2) {
                    def detected_lane_r2 = extractLaneFromFilename(fastq_2.toString())
                    if (detected_lane != detected_lane_r2) {
                        log.error("""
                            Lane mismatch for paired reads of sample '${meta.id}'
                            Paired-end reads must originate from the same sequencing run.

                            Read 1: ${fastq_1}
                            Lane: ${detected_lane}

                            Read 2: ${fastq_2}
                            Lane: ${detected_lane_r2}
                            """.stripIndent())
                        error("Please check your samplesheet for mismatched FASTQ file pairs.")
                    }
                }

                meta = meta + [lane: detected_lane]
            }

            // if lane is manually provided, validate that it matches FASTQ filename (optional strict mode)
            else if ( (params.strict_lane_detection || params.auto_detect_lanes) && meta.lane) {
                def detected_lane = extractLaneFromFilename(fastq_1.toString())
                if (detected_lane && meta.lane != detected_lane) {
                    log.warn("""
                        Samplesheet provided lane '${meta.lane}'
                        does not match detected lane '${detected_lane}'
                        extracted from file name '${fastq_1}'.
                        Using provided value in samplesheet: '${meta.lane}'.
                    """.stripIndent())
                }

                // check if detectes lanes match for paired reads
                if (fastq_2) {
                    def detected_lane_r2 = extractLaneFromFilename(fastq_2.toString())
                    if (detected_lane &&
                        detected_lane_r2 &&
                        detected_lane != detected_lane_r2) {
                        log.error("""
                            Lane mismatch for paired reads of sample '${meta.id}'
                            Paired-end reads must originate from the same sequencing run.

                            Provided lane: ${meta.lane}

                            Read 1: ${fastq_1}
                            Detected, but ignored lane: ${detected_lane}

                            Read 2: ${fastq_2}
                            Detected lane: ${detected_lane_r2}
                            """.stripIndent())
                        error("Please check your samplesheet for mismatched FASTQ file pairs.")
                    }
                }

            }
            // auto-detect flowcell if not provided - used for checking uniqueness and read group construction
            if (params.auto_detect_flowcells && !meta.flowcell) {
                def detected_flowcell = extractFlowcellFromFastq(fastq_1)
                if (detected_flowcell) {
                    log.warn("Auto-detected flowcell ${detected_flowcell} for sample ${meta.id} with fastq file(s) ${fastq_1} ${fastq_2}")
                }
                else {
                    log.warn("Proceeding without flowcell info for ${meta.id}, but this could lead to duplicate input warnings or incorrect results.")
                }
                // check if flowcells match for paired reads
                if (fastq_2) {
                    def detected_flowcell_r2 = extractFlowcellFromFastq(fastq_2)
                    if (detected_flowcell != detected_flowcell_r2) {
                        log.error("""
                            Flowcell ID mismatch for paired reads of sample '${meta.id}'
                            Paired-end reads must originate from the same sequencing run.

                            Read 1: ${fastq_1}
                            Flowcell: ${detected_flowcell}

                            Read 2: ${fastq_2}
                            Flowcell: ${detected_flowcell_r2}
                            """.stripIndent())
                        error("Please check your samplesheet for mismatched FASTQ file pairs.")
                    }
                }
                meta = meta + [flowcell: detected_flowcell]
            }
            // If flowcell is manually provided, validate it matches FASTQ headers (optional strict mode)
            else if (params.strict_flowcell_detection && params.auto_detect_flowcells && meta.flowcell) {
                def detected_flowcell = extractFlowcellFromFastq(fastq_1)
                if (detected_flowcell && meta.flowcell != detected_flowcell) {
                    log.warn("""
                        Samplesheet provided flowcell '${meta.flowcell}'
                        does not match detected flowcell '${detected_flowcell}'
                        extracted from FASTQ header in '${fastq_1}'.
                        Using provided value in samplesheet: '${meta.flowcell}'.
                    """.stripIndent())
                }
                if (fastq_2) {
                    def detected_flowcell_r2 = extractFlowcellFromFastq(fastq_2)
                    if (detected_flowcell &&
                        detected_flowcell_r2 &&
                        detected_flowcell != detected_flowcell_r2) {
                        log.error("""
                            Flowcell ID mismatch for paired reads of sample '${meta.id}'
                            Paired-end reads must originate from the same sequencing run.

                            Provided flowcell: ${meta.flowcell}

                            Read 1: ${fastq_1}
                            Detected, but ignored flowcell: ${detected_flowcell}

                            Read 2: ${fastq_2}
                            Detected flowcell: ${detected_flowcell_r2}
                            """.stripIndent())
                        error("Please check your samplesheet for mismatched FASTQ file pairs.")
                    }
                }
            }
            // handle single and paired-end fastq files
            if (!fastq_2) {
                return [ meta.id, [ meta + [ single_end:true, sample: "${meta.id}" ], [ fastq_1 ] ] ]
            } else {
                return [ meta.id, [ meta + [ single_end:false, sample: "${meta.id}" ], [ fastq_1, fastq_2 ] ] ]
            }
            // TODO compare with nf-core template approach, which uses fewer nested lists + flatten.
        }
        .tap{ ch_samples } // save sample-wise channel for later re-use
        // [sample_1, [[id:sample_1, lane:1, library:null, flowcell:232KF7LT4, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]]
        // group by sample id
        .groupTuple()
        // [sample_1, [[[id:sample_1, lane:1, library:null, flowcell:232KF7LT4, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]], [[id:sample_1, lane:2, library:null, flowcell:232WTNLT4, single_end:false, sample:sample_1], [/path/to/sample_1_L002_R1_001.fastq.gz, /path/to/sample_1_L002_R2_001.fastq.gz]]]]
        .map { validateInputSamplesheet(it) }
        // calculate number of fastq (pairs) per sample
        .map { sample, ch_items ->
            def metas = ch_items.collect { it[0] }
            def num_entries = ch_items.size()
            def has_multiple_lanes = metas.collect { it.lane }.findAll { it != null }.unique().size() > 1
            def has_multiple_libraries = metas.collect { it.library }.findAll { it != null }.unique().size() > 1
            def has_multiple_flowcells = metas.collect { it.flowcell }.findAll { it != null }.unique().size() > 1
            return [ sample, num_entries, has_multiple_lanes, has_multiple_libraries, has_multiple_flowcells ]
        }
        // [sample_1, 2, true, false, true]
        .combine(ch_samples, by: 0)
        // [sample_1, 2, true, false, true, [[id:sample_1, lane:1, library:null, flowcell:232KF7LT4, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]]
        // [sample_1, 2, true, false, true, [[id:sample_1, lane:2, library:null, flowcell:232WTNLT4, single_end:false, sample:sample_1], [/path/to/sample_1_L002_R1_001.fastq.gz, /path/to/sample_1_L002_R2_001.fastq.gz]]]
        .map { _sample, num_entries, multi_lane, multi_lib, multi_fc, ch_item ->
            def ( meta, fastqs ) = ch_item

            // build unique identifier based on sample name, library, lane and flowcell
            def id_elements = [sanitizeMetadataValue(meta.sample)]
            if (meta.library)
                id_elements.add("LIB_${sanitizeMetadataValue(meta.library)}")
            if (meta.lane)
                id_elements.add("LANE_${sanitizeMetadataValue(meta.lane)}")
            if (meta.flowcell)
                id_elements.add("FC_${sanitizeMetadataValue(meta.flowcell)}")
            meta = meta + [
                id: id_elements.join('-'),
                num_entries: num_entries.toInteger(),
                multiple_lanes: multi_lane,
                multiple_libraries: multi_lib,
                multiple_flowcells: multi_fc
            ]

            return [ meta, fastqs ]
        }
        // [[id:sample_1-LANE_1-FC_232KF7LT4, lane:1, library:null, flowcell:232KF7LT4, single_end:false, sample:sample_1, num_entries:2, multiple_lanes:true, multiple_libraries:false], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]]
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
    monochrome_logs // boolean: Disable ANSI colour codes in log output

    main:

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)

    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs for common issues: https://nf-co.re/docs/running/troubleshooting"
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
    // [sample_1, [[[id:sample_1, lane:1, library:null, flowcell:232KF7LT4, single_end:false, sample:sample_1], [/path/to/sample_1_L001_R1_001.fastq.gz, /path/to/sample_1_L001_R2_001.fastq.gz]], [[id:sample_1, lane:2, library:null, flowcell:232WTNLT4, single_end:false, sample:sample_1], [/path/to/sample_1_L002_R1_001.fastq.gz, /path/to/sample_1_L002_R2_001.fastq.gz]]]]
    def (sample, items) = input
    def metas = items.collect { it[0] }

    // Only validate if there are multiple fastq (pair) entries for the same sample
    if (metas.size > 1) {
        // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
        def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
        if (!endedness_ok) {
            log.error("Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
            error("Please check your samplesheet for consistent single/paired-end entries per sample.")
        }

        // Check fi samples can be distinguished based on lane/flowcell/library information
        def identifiers = metas.collect { "${it.flowcell ?: 'NA'}_${it.lane ?: 'NA'}_${it.library ?: 'NA'}" }
        if (identifiers.size() != identifiers.unique().size()) {
            log.error("""
                Sample '${sample}' has duplicate flowcell/lane/library combinations.
                Each sample must have a unique combination of flowcell, lane and library identifiers.
                Flowcell IDs were automatically extracted from FASTQ headers.
            """.stripIndent())
            error("Please check samplesheet for duplicate entries.")
        }

        // Could in theory be replaced by adding the following to schema_input.json:
        // "uniqueEntries": ["sample", "lane", "library"]
        // Should be placed just after the items declaration (i.e. level 0).
        // However, error message will be less detailed.
        // Also needs verification whether shared null/[] entries would trigger an error or not.
        // This also won't take into account flowcell IDs.

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
        } else {
            def msg = "Failed to extract lane info from FASTQ file name: ${filename}"
            if (params.strict_lane_detection) {
                error(msg)
            } else {
                log.warn(msg)
            }
        }
        return null // return nothing if no pattern is found
    }
}

// Parse first line of a FASTQ file, return the flowcell ID
// Adapted from nf-core/sarek https://github.com/nf-core/sarek/blob/5cc30494a6b8e7e53be64d308b582190ca7d2585/workflows/sarek/main.nf#L953
// Was originally done at run-time to construct read groups, but was moved
// to input validation for improved uniqueness validation.
// Now it is possible to use identical sample names and lanes for runs on different flowcells
def extractFlowcellFromFastq(path) {
    // expected formats (https://en.wikipedia.org/wiki/FASTQ_format#Illumina_sequence_identifiers)
    //
    // CASAVA 1.8 - seven fields
    // @A00379:443:HJMTYDSX2:4:1101:4029:1078 1:N:0:TCGTGGAT+ACCTAGAC
    // @<instrument>:<run-id>:<flowcell-id>:<lane>:<tile>:<x-coord>:<y-coord> <read-pair>:<filtered>:<control-bits>:<index-sequence>
    //
    // Older Illumina format - five fields
    // @HWUSI-EAS100R:6:73:941:1973#0/1
    // @<machine_id>:<lane>:<tile>:<x_coord>:<y_coord>#<index>/<read>
    //
    // In SRA/ENI fastq files, the identifier often starts with an accession, followed by a space, followed by the actual fastq identifier
    // For the five fields format, this can lead to a flowcell id containing a space

    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)

    String fcid = null
    String line

    try {
        path.withInputStream {
            InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
            Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
            BufferedReader buffered = new BufferedReader(decoder)
            line = buffered.readLine()
        }
    } catch (Exception e) {
        log.error("ERROR: Could not extract flowcell ID from ${path}: ${e.message}")
        error("Please verify that the file is a valid gzipped FASTQ file.")

    }

    if (!line || !line.startsWith('@')) {   // ~ !line?.startsWith('@')
        log.error("ERROR: Invalid FASTQ header in ${path}: ${line}")
        error("Please verify that this is a valid FASTQ file with headers that start with '@'.")
    }

    // remove leading @ and strip whitespace
    line = line.substring(1).trim()

    // split on whitespace to account for SRA identifiers
    def tokens = line.split(/\s+/)
    def header = tokens.find { it.count(':') >= 4 }

    if (!header) {
        def msg = """
            Could not find Illumina-style header in FASTQ ${path}:
            ${line}
        """.stripIndent()

        if (params.strict_flowcell_detection) {
            error(msg)
        } else {
            log.warn(msg)
            return null
        }
    }

    def fields = header.split(':')

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format
        fcid = fields[2]
    } else if (fields.size() == 5) {
        // Old Illumina format
        fcid = fields[0]
    }

    if (!fcid) {
        def msg = """
            Could not determine flowcell ID from FASTQ header ${path}:
            ${line}

            Expected one of:
              - CASAVA 1.8+ format
              - Legacy Illumina format
              - SRA-reformatted Illumina format

            Flowcell auto-detection will be skipped.
        """.stripIndent()
        if (params.strict_flowcell_detection) {
            error(msg)
        } else {
            log.warn(msg)
            return null
        }
    }

    return fcid
}

def sanitizeMetadataValue(value) {
    return value?.toString()?.replaceAll(/[^A-Za-z0-9._-]/, '_')
}

// TODO: optionally just return warnings instead of stopping the pipeline. Non-unique entries should still be caught further downstream and this would allow non-standard files to be processed if there is just a single run per sample (or even just no repeated lanes/lib combinations).

// Test cases
// assert extractLaneFromFilename('sample_L001_R1.fastq.gz') == '001'
// assert extractLaneFromFilename('sample_L1_R1.fastq.gz') == '001'
// assert extractLaneFromFilename('sample_lane2_R1.fastq.gz') == '002'
// assert extractLaneFromFilename('sample_100_R1.fastq.gz') == null  // Not a lane!

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "fastp (Chen et al. 2018),",
            "FastQ Screen (Wingett and Andrews 2018),",
            "deacon (Constantinides et al. 2025),",
            "bbsplit (Bushnell, n.d.),",
            "SeqKit2 (Shen et al. 2024),",
            "BWA-MEM (Li 2013),",
            "Picard Toolkit (Broad Institute 2019),",
            "SAMtools (Li et al. 2009),",
            "Mosdepth (Pedersen and Quinlan 2018),",
            "GATK (Van der Auwera and O'Connor 2020; Poplin et al. 2017),",
            "SnpEff (Cingolani et al. 2012),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews, Simon. 2010. FastQC: A Quality Control Tool for High Throughput Sequence Data. Released. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.</li>",
            "<li>Chen, Shifu, Yanqing Zhou, Yaru Chen, and Jia Gu. 2018. “Fastp: An Ultra-Fast All-in-One FASTQ Preprocessor.” Bioinformatics 34 (17): i884–90. https://doi.org/10.1093/bioinformatics/bty560.</li>",
            "<li>Wingett, Steven W., and Simon Andrews. 2018. “FastQ Screen: A Tool for Multi-Genome Mapping and Quality Control.” F1000Research 7 (September): 1338. https://doi.org/10.12688/f1000research.15931.2.</li>",
            "<li>Constantinides, Bede, John Lees, and Derrick W. Crook. 2025. “Deacon: Fast Sequence Filtering and Contaminant Depletion.” Preprint, Bioinformatics, June 12. https://doi.org/10.1101/2025.06.09.658732.</li>",
            "<li>Bushnell, Brian. n.d. Bbmap/Bbsplit. https://sourceforge.net/projects/bbmap/.</li>",
            "<li>Shen, Wei, Botond Sipos, and Liuyang Zhao. 2024. “SeqKit2: A Swiss Army Knife for Sequence and Alignment Processing.” iMeta 3 (3): e191. https://doi.org/10.1002/imt2.191.</li>",
            "<li>Li, Heng. 2013. “Aligning Sequence Reads, Clone Sequences and Assembly Contigs with BWA-MEM.” arXiv:1303.3997. Version 2. Preprint, arXiv. https://doi.org/10.48550/ARXIV.1303.3997.</li>",
            "<li>Broad Institute. 2019. Picard Toolkit. Released. https://broadinstitute.github.io/picard/.</li>",
            "<li>Li, Heng, Bob Handsaker, Alec Wysoker, et al. 2009. “The Sequence Alignment/Map Format and SAMtools.” Bioinformatics 25 (16): 2078–79. https://doi.org/10.1093/bioinformatics/btp352.</li>",
            "<li>Pedersen, Brent S., and Aaron R. Quinlan. 2018. “Mosdepth: Quick Coverage Calculation for Genomes and Exomes.” Bioinformatics 34 (5): 867–68. https://doi.org/10.1093/bioinformatics/btx699.</li>",
            "<li>Van Der Auwera, Geraldine A., and Brian D. O’Connor. 2020. Genomics in the Cloud: Using Docker, GATK, and WDL in Terra (1st Edition). 1st ed. O’Reilly Media. https://www.oreilly.com/library/view/genomics-in-the/9781491975183/.</li>",
            "<li>Poplin, Ryan, Valentin Ruano-Rubio, Mark A. DePristo, et al. 2017. “Scaling Accurate Genetic Variant Discovery to Tens of Thousands of Samples.” Preprint. Genomics, November 14. https://doi.org/10.1101/201178.</li>",
            "<li>Cingolani, Pablo, Adrian Platts, Le Lily Wang, et al. 2012. “A Program for Annotating and Predicting the Effects of Single Nucleotide Polymorphisms, SnpEff: SNPs in the Genome of Drosophila Melanogaster Strain w1118 ; Iso-2; Iso-3.” Fly 6 (2): 80–92. https://doi.org/10.4161/fly.19695.</li>",
            "<li>Ewels, Philip, Måns Magnusson, Sverker Lundin, and Max Käller. 2016. “MultiQC: Summarize Analysis Results for Multiple Tools and Samples in a Single Report.” Bioinformatics 32 (19): 3047–48. https://doi.org/10.1093/bioinformatics/btw354.</li>"
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

    meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
