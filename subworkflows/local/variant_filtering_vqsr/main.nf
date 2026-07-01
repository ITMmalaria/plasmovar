//
// SUBWORKFLOW: VARIANT_FILTERING_VQSR
//

// For more info on order of steps, see:
//  https://github.com/broadgsa/gatk/blob/master/doc_archive/tutorials/(howto)_Recalibrate_variant_quality_scores_=_run_VQSR.md,
//  https://github.com/broadgsa/gatk/blob/master/doc_archive/faqs/,Which_training_sets___arguments_should_I_use_for_running_VQSR%3F.md, and
//  https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering,
// or the MalariaGEN Pf8 pipeline: https://github.com/malariagen/malariagen-pf8-snp-indel-calling/blob/529fe6b59bbf14a99d8c46264f6a38c4761a0ffa/main.nf#L134

include { GATK4_MERGEVCFS as GATHER_VCFS_BEFORE_VQSR                   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MAKESITESONLYVCF                                       } from '../../../modules/local/gatk4/makesitesonlyvcf/main'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_SNP   } from '../../../modules/local/gatk4/variantrecalibrator'
include { GATK4_VARIANTRECALIBRATOR as GATK4_VARIANTRECALIBRATOR_INDEL } from '../../../modules/local/gatk4/variantrecalibrator'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_SNP                       } from '../../../modules/nf-core/gatk4/applyvqsr/main'
include { GATK4_APPLYVQSR as GATK4_APPLYVQSR_INDEL                     } from '../../../modules/nf-core/gatk4/applyvqsr/main'

workflow VARIANT_FILTERING_VQSR {

    take:
    ch_vcf_by_interval  // [[meta], vcf.gz, vcf.gz.tbi] - 1 per interval
    ch_ref_fasta        // [[meta], fasta]
    ch_ref_fai          // [[meta], fasta.fai]
    ch_ref_dict         // [[meta], dict]
    // snp_args_file
    // indel_args_file

    main:

    // Prepare for VCF gathering: collect all interval VCFs without index
    // Group by genome_id to merge all intervals for that genome
    // Sort on interval index for deterministic channel order
    ch_vcfs_to_gather = ch_vcf_by_interval                              // [ meta, vcf, tbi ]                                       (1 per interval)
        .map { meta, vcf, _tbi ->
            [ meta.genome_id, [index: meta.interval_index, vcf: vcf] ]  // [ meta, [ index, vcf ] ]                                 (1 per interval)
        }
        .groupTuple(by: 0)                                              // [ meta.genome_id, [ [index, vcf], [index_vcf], ... ] ]   (1 element with all per-interval files)
        .map { genome_id, vcfs ->
            // sort on interval_list index to make channel inputs for gathering interval VCFs more deterministic
            def sorted_vcfs = vcfs
                .sort { it.index }
                .collect { it.vcf }
            def final_meta = [ id: "${genome_id}_merged", genome_id: genome_id ]
            [ final_meta, sorted_vcfs ]                                 // [ meta, [ vcfs ] ]                                       (1 element with all per-interval files)
        }

    // Gather interval VCFs into single cohort VCF
    // Order is handled automatically by the sequence dictionary
    // (unlike bcftools, see https://nf-co.re/subworkflows/vcf_gather_bcftools)
    GATHER_VCFS_BEFORE_VQSR(
        ch_vcfs_to_gather,
        ch_ref_dict
    )
    ch_gathered_vcf_tbi = GATHER_VCFS_BEFORE_VQSR.out.vcf.join(GATHER_VCFS_BEFORE_VQSR.out.tbi)

    // TODO: optional add ExcessHet filter for large cohorts. See https://gatk.broadinstitute.org/hc/en-us/articles/360035531112--How-to-Filter-variants-either-with-VQSR-or-by-hard-filtering

    // Create sites-only VCF - VQSR modelling does not require sample-level annotations
    GATK4_MAKESITESONLYVCF(ch_gathered_vcf_tbi)
    ch_gathered_sites_only_vcf_tbi_for_vqsr_modelling = GATK4_MAKESITESONLYVCF.out.vcf.join(GATK4_MAKESITESONLYVCF.out.tbi)

    // Create resource channels
    ch_snp_vcfs = channel.fromPath(params.vqsr_snp_resource_vcfs, checkIfExists: true)
    ch_snp_tbis = channel.fromPath(params.vqsr_snp_resource_vcfs.collect { tbi -> "${tbi}.tbi" }, checkIfExists: true)

    ch_indel_vcfs = channel.fromPath(params.vqsr_indel_resource_vcfs, checkIfExists: true)
    ch_indel_tbis = channel.fromPath(params.vqsr_indel_resource_vcfs.collect { tbi -> "${tbi}.tbi" }, checkIfExists: true)

    // Build VQSR model for SNPs - uses sites-only vcf input
    GATK4_VARIANTRECALIBRATOR_SNP(
        ch_gathered_sites_only_vcf_tbi_for_vqsr_modelling,
        'SNP',
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict },
        ch_snp_vcfs.collect(),
        ch_snp_tbis.collect(),
        params.vqsr_snp_resource_labels
    )

    // Build VQSR model for indels - uses sites-only VCF input
    GATK4_VARIANTRECALIBRATOR_INDEL(
        ch_gathered_sites_only_vcf_tbi_for_vqsr_modelling,
        'INDEL',
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict },
        ch_indel_vcfs.collect(),
        ch_indel_tbis.collect(),
        params.vqsr_indel_resource_labels
    )

    // Apply VQSR to SNPs - use original VCF with sample-level genotype annotations
    GATK4_APPLYVQSR_SNP(
        ch_gathered_vcf_tbi
            .combine(GATK4_VARIANTRECALIBRATOR_SNP.out.recal.map { _meta, recal -> recal })
            .combine(GATK4_VARIANTRECALIBRATOR_SNP.out.idx.map { _meta, idx -> idx })
            .combine(GATK4_VARIANTRECALIBRATOR_SNP.out.tranches.map { _meta, tranches -> tranches }),
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict }
    )
    ch_gathered_vcf_tbi_snp_recal = GATK4_APPLYVQSR_SNP.out.vcf.join(GATK4_APPLYVQSR_SNP.out.tbi)

    // Apply VQSR to INDELs - uses recalibrated SNP VCF as input
    GATK4_APPLYVQSR_INDEL(
        ch_gathered_vcf_tbi_snp_recal
            .combine(GATK4_VARIANTRECALIBRATOR_INDEL.out.recal.map { _meta, recal -> recal })
            .combine(GATK4_VARIANTRECALIBRATOR_INDEL.out.idx.map { _meta, idx -> idx })
            .combine(GATK4_VARIANTRECALIBRATOR_INDEL.out.tranches.map { _meta, tranches -> tranches }),
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict }
    )

    emit:
    vcf_filter_added     = GATK4_APPLYVQSR_INDEL.out.vcf    // [[meta], sorted.vcf.gz] - with filter tags/annotations, no omitted sites
    vcf_filter_added_tbi = GATK4_APPLYVQSR_INDEL.out.tbi    // [[meta], sorted.vcf.gz.tbi]
}
