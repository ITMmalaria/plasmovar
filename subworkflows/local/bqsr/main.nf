//
// SUBWORKFLOW: Base Quality Score Recalibration (BQSR)
//
// Utilises scatter-gather approach, i.e. input should be per interval channel
// Outputs calibrated bam files for each sample, scattered per interval
// (note that bam file is not physically split up, but each element holds a specific interval_list file)
//

// TODO: add fall-back to uncalibrated bam in case of applyBQSR errors

include { GATK4_BEDTOINTERVALLIST } from '../../../modules/nf-core/gatk4/bedtointervallist/main'
include { GATK4_INTERVALLISTTOOLS } from '../../../modules/nf-core/gatk4/intervallisttools/main'
include { GATK4_BASERECALIBRATOR  } from '../../../modules/nf-core/gatk4/baserecalibrator/main'
include { GATK4_GATHERBQSRREPORTS } from '../../../modules/nf-core/gatk4/gatherbqsrreports/main'
include { GATK4_APPLYBQSR         } from '../../../modules/nf-core/gatk4/applybqsr/main'

workflow BQSR {

    take:
    ch_bam_bai_intervals    // [ [combined_meta], bam, bai, interval ]
    ch_ref_fasta            // [ [meta], fasta ]
    ch_ref_fai              // [ [meta], fai ]
    ch_ref_dict             // [ [meta], dict]
    known_sites_vcf         // path to vcf file
    known_sites_tbi         // path to vcf tbi file
    scatter                 // boolean
    bed                     // [bed]

    main:

    // Prepare channels with vcf and tbi files of known sites
    vcf = file(known_sites_vcf, checkIfExists: true)
    ch_known_sites_vcf = channel.value([[id: "known_sites"], vcf])
    tbi = file(known_sites_tbi, checkIfExists: true)
    ch_known_sites_tbi = channel.value([[id: "known_sites"], tbi])

    // Gather bam files to sample-level
    // This channel will be used by ApplyBQSR (which should never run in scattered mode)
    // and also provide the input for BaseRecalibrator when running in non-scatter mode
    // The latter can help avoid problems if known site vcf contains fewer regions than
    // in the reference genome. Likewise for the custom bqsr_bed mode.
    ch_bam_gathered = ch_bam_bai_intervals
        .map{ meta, bam, bai, _interval ->
            def sample_meta = (meta - [ interval_name: meta.interval_name ]) + [ id: meta.sample ]
            tuple(meta.sample, sample_meta, bam, bai)
        }
        // reduce to a single bam per sample now that intervals are removed
        .distinct{ it[0] }  // distinct by sample name
        .map { _sample, sample_meta, bam, bai -> tuple(sample_meta, bam, bai) }

    // Prepare channel for BaseRecalibrator input
    // In non-scatter mode, simply use the sample-level bam channel
    if (!scatter) {
        // add empty element for GATK4_BASERECALIBRATOR
        ch_bam_baserecalibrator = ch_bam_gathered.map { meta, bam, bai -> [meta, bam, bai, []]}
    }
    // If a custom interval (bed) file is provided, convert it to an interval_list and re-create interval-scattered bam channel
    // TODO: these steps could be combined into a subworkflow, which could then be used both here and for the main GATK interval creation step
    // TODO: import modules as *_BQSR, so that separate modules.config settings can be used for e.g. publishDir. Or ensure different meta.id.
    else if (bed) {
        // Convert BED to GATK IntervalList
        bed_file = file(bed, checkIfExists: true)
        ch_ref_bed = channel.value([[id: bed_file.simpleName], bed_file])
        GATK4_BEDTOINTERVALLIST(ch_ref_bed, ch_ref_dict)

        // Combine intervals into limited number of separate interval_list files for scatter-gather parallel processing
        GATK4_INTERVALLISTTOOLS(GATK4_BEDTOINTERVALLIST.out.interval_list)
        ch_intervals = GATK4_INTERVALLISTTOOLS.out.interval_list
            // [ [interval_genome_meta.id], [interval_list, interval_list, ...] ] (single element)
            .transpose()
            // [ [interval_genome_meta.id], interval_list ] (multiple elements)

        // Combine with custom bqsr interval lists
        ch_bam_baserecalibrator = ch_bam_gathered
            .combine(ch_intervals)
            .map { bam_meta, bam, bai, interval_genome_meta, interval ->
                def combined_meta = bam_meta + [
                    // add unique id per sample/interval combination
                    id: "${bam_meta.id}_${interval.simpleName}",
                    genome_id: interval_genome_meta.id,
                    interval_name: interval.simpleName
                ]
                [ combined_meta, bam, bai, interval ]
            }
            // [ [combined_meta.id, combined_meta.sample, combined_meta.num_entries, combined_meta.multiple_lanes, combined_meta.multiple_libraries, combined_meta.data_type, combined_meta.genome_id, combined_meta.interval_name], bam, bai, interval_list ]
    }
    // In normal scatter mode, use the original gatk intervals
    else {
        ch_bam_baserecalibrator = ch_bam_bai_intervals
    }

    // Run BaseRecalibrator
    GATK4_BASERECALIBRATOR(
        ch_bam_baserecalibrator,
        ch_ref_fasta,
        ch_ref_fai,
        ch_ref_dict,
        ch_known_sites_vcf,
        ch_known_sites_tbi
    )

    // Gather scattered tables into a single table for ApplyBQSR when using scatter mode (either custom or original gatk intervals)
    if (scatter || bed) {
        ch_table = GATK4_BASERECALIBRATOR.out.table  // [[meta], table]
            // 1) remove interval-specific meta field
            // 2) re-create original sample id without interval suffix
            // to make meta identical for all tables belonging to a particular sample
            // and allow grouping of table file paths into a single element (per sample)
            .map { meta, table ->
                def new_meta = (meta - [ interval_name: meta.interval_name ]) + [ id: meta.sample ]
                // If meta were to change in the future, it would be more robust to only retain specific meta fields instead of omitting the ones that might interfere
                // new_meta = meta.subMap(['id', 'sample', ...]) + [id: meta.sample]
                tuple(new_meta, table)
            }
            .groupTuple()

        GATK4_GATHERBQSRREPORTS(ch_table)

        // Combine gathered table with interval-level bam channel
        // Key by sample id explicitly, because full meta differs (interval_name is missing from table channel)
        ch_table = GATK4_GATHERBQSRREPORTS.out.table.map { meta, table -> [meta.sample, table] }
    }
    // In non-scatter mode, combine the single recalibration table with sample-level bam channel
    else if (!scatter) {
        ch_table = GATK4_BASERECALIBRATOR.out.table.map { meta, table -> [meta.sample, table] }
    }

    // Combine sample-level bam channel with the single recalibration table
    // Alternatively, a join (not combine!) on interval-level bam channel could be used too (join reduces to a single unique sample)
    // Note that the metamap of ch_bam_gathered contains sample-level info, e.g.
    // id = meta.sample, instead of meta.sample-interval and there is no meta.interval_name
    ch_bam_bai_table = ch_bam_gathered
        // key by sample name for joining with recalibration table
        .map { sample_meta, bam, bai -> [sample_meta.sample, sample_meta, bam, bai] }
        .join(ch_table, by: 0)
        // Pass empty list for the expected intervals input of GATK4_APPLYBQSR,
        // since it should not be run in scatter mode according to GATK.
        // See comments at bottom of file for more info.
        .map { _sample_key, sample_meta, bam, bai, table -> tuple(sample_meta, bam, bai, table, []) }

    // Apply recalibration to per-sample BAM files without scattering
    GATK4_APPLYBQSR(
        ch_bam_bai_table,
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict }
    )
    // !NOTE: applybqsr might fail with read group not found errors if there are too few reads mapped/too low quality
    // org.broadinstitute.hellbender.exceptions.GATKException: Read group pf1.null not found in the recalibration table.Set the allow-missing-read-group command line argument to ignore this error.
    // see: https://gatk.broadinstitute.org/hc/en-us/community/posts/360076232572--Repost-ReadGroup-missing-when-Applying-BQSR
    // also observed in pf8 pipeline: https://github.com/malariagen/malariagen-pf8-snp-indel-calling/blob/529fe6b59bbf14a99d8c46264f6a38c4761a0ffa/modules/bqsr.nf#L50

    // TODO: add AnalyzeCovariates + multi-pass bqsr
    // cf. https://github.com/nf-core/genomicrelatedness/blob/dev/subworkflows/local/base_quality_score_recalibration/main.nf
    // https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/
    // > BQSR is performed twice. The second pass is optional, only required to produce a recalibration report.

    emit:
    bam_recalibrated = GATK4_APPLYBQSR.out.bam  // [[meta], bam]
    bai_recalibrated = GATK4_APPLYBQSR.out.bai  // [[meta], bai]
}

// Info on scatter-gather approach for BQSR

// https://gatk.broadinstitute.org/hc/en-us/community/posts/360077410652-What-interval-should-be-used-when-doing-the-BaseRecalibrator-for-exome-data

// I contacted the workflow maintainers and here is their response to the question why in the workflow BaseRecalibrator is given automatically generated interval files from CreateSequenceGroupingTSV task using the reference dictionary instead of the intervals specific to their exome.

//     It’s important to note that the scattered BQSR reports are gathered and combined before ApplyBQSR is run, so doing this in the scattered approach, where our scattering intervals are whole chromosomes, gives exactly the same results as if it were run unscattered.  In that sense, assuming the scattering intervals fully tile the entire reference, the particular choice of the scattering intervals has no effect on the output, and is selected based on resource and time considerations.

// Furthermore regarding whether BaseRecalibrator should be run subset to only the regions around the targets, regardless of whether it is being scatter-gathered or not.

//     This is a more complex question.  The effect of subsetting to the target regions for BaseRecalibrator would be that the base error model learned by BaseRecalibrator would only be based on reads in/near the target region.  This could be a reasonable choice, since reads aligned outside of these regions would be more likely to be mapping errors, and so could artificially inflate the error rate learned by BaseRecalibrator.  In practicality, however, I would expect the difference from subsetting to the target regions to be negligible in most data.  I would only expect a noticeable difference in data with a large fraction of off target reads, and such data would be flagged as problematic by hybrid selection related QC metrics.  So in the end, I would say the advantage of maintaining a single code path used for both exomes and whole genomes for BaseRecalibrator carries the day, since the potential risks are very minimal.

//     Also note, ApplyBQSR must ALWAYS be run over the full genome + unmapped (whether in a single job, or a scatter-gather).  Otherwise, some reads can lose their mates in the output data.

// https://github.com/nf-core/sarek/issues/445
// By definition, exome sequencing and other targeted sequencing data don’t cover the entire genome, so most analyses can be
// restricted to just the capture targets (genes or exons) to save processing time and enable scatter gather parallelism. In
// addition, there are some processing steps, such as BQSR, that should be restricted to the capture targets in order to
// eliminate off-target sequencing data, which is uninformative and is a source of noise.

// https://gatk.broadinstitute.org/hc/en-us/articles/360035889551-When-should-I-restrict-my-analysis-to-specific-intervals
// https://gatk.broadinstitute.org/hc/en-us/community/posts/4403947098395-What-intervals-to-use-in-HaplotypeCaller-GenotypeGVCF-for-exome-data-sequenced-with-different-capture-kits

// https://github.com/nf-core/sarek/issues/1772
// https://gatk.broadinstitute.org/hc/en-us/community/posts/360077103652-output-BAM-form-ApplyBQSR-query-using-interval-option
// You should not use intervals with ApplyBQSR because you want to recalibrate all of the reads. You can introduce artifacts if you run with intervals in ApplyBQSR. There is an in depth discussion on the forum at this link: https://gatk.broadinstitute.org/hc/en-us/community/posts/360074603551-ValidateSamFile-New-error-after-running-BQSR

// In the ApplyBQSR docs, the -L argument is not a recommended parameter.

// https://gatk.broadinstitute.org/hc/en-us/community/posts/360074603551-ValidateSamFile-New-error-after-running-BQSR
// With BQSR, intervals can be used while running BaseRecalibrator but not ApplyBQSR. If you use intervals, you will need to add an extra step after BaseRecalibrator to combine all the recalibration tables. ApplyBQSR is meant to be run on all the reads that contribute to the model so it cannot be run using intervals.
