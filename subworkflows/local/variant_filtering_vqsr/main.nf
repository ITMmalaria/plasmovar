//
// SUBWORKFLOW: VARIANT_FILTERING_VQSR
//

include { GATK4_MERGEVCFS as GATHER_VCFS_BEFORE_VQSR                   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
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
    ch_vcfs_to_gather = ch_vcf_by_interval
        // [ meta, vcf ]    (1 per interval)
        .map { meta, vcf, _tbi ->
            [ meta.genome_id, vcf ]
        }
        .groupTuple(by: 0)
        // [ meta.genome_id, [vcfs] ] (1 element with all interval files)
        .map { genome_id, vcfs ->
            def final_meta = [ id: "${genome_id}_merged", genome_id: genome_id ]
            [ final_meta, vcfs ]
        }

    // Gather interval VCFs into single cohort VCF
    // Order is handled automatically by the sequence dictionary
    // (unlike bcftools, see https://nf-co.re/subworkflows/vcf_gather_bcftools)
    GATHER_VCFS_BEFORE_VQSR(
        ch_vcfs_to_gather,
        ch_ref_dict
    )
    ch_gathered_vcf_tbi_for_vqsr = GATHER_VCFS_BEFORE_VQSR.out.vcf.join(GATHER_VCFS_BEFORE_VQSR.out.tbi)

    // Create resource channels
    ch_snp_vcfs = channel.fromPath(params.vqsr_snp_resource_vcfs, checkIfExists: true)
    ch_snp_tbis = channel.fromPath(params.vqsr_snp_resource_vcfs.collect { tbi -> "${tbi}.tbi" }, checkIfExists: true)

    ch_indel_vcfs = channel.fromPath(params.vqsr_indel_resource_vcfs, checkIfExists: true)
    ch_indel_tbis = channel.fromPath(params.vqsr_indel_resource_vcfs.collect { tbi -> "${tbi}.tbi" }, checkIfExists: true)

    // Build VQSR model for SNP
    GATK4_VARIANTRECALIBRATOR_SNP(
        ch_gathered_vcf_tbi_for_vqsr,
        'SNP',
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict },
        ch_snp_vcfs.collect(),
        ch_snp_tbis.collect(),
        params.vqsr_snp_resource_labels
    )

    // Apply VQSR for SNPs
    GATK4_APPLYVQSR_SNP(
        ch_gathered_vcf_tbi_for_vqsr
            .combine(GATK4_VARIANTRECALIBRATOR_SNP.out.recal.map { _meta, recal -> recal })
            .combine(GATK4_VARIANTRECALIBRATOR_SNP.out.idx.map { _meta, idx -> idx })
            .combine(GATK4_VARIANTRECALIBRATOR_SNP.out.tranches.map { _meta, tranches -> tranches }),
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict }
    )

    // Build VQSR model for indel using the recalibrated SNP VCF
    // For more info on order of steps, see https://github.com/broadgsa/gatk/blob/master/doc_archive/tutorials/(howto)_Recalibrate_variant_quality_scores_=_run_VQSR.md
    GATK4_VARIANTRECALIBRATOR_INDEL(
        GATK4_APPLYVQSR_SNP.out.vcf.join(GATK4_APPLYVQSR_SNP.out.tbi),
        'INDEL',
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict },
        ch_indel_vcfs.collect(),
        ch_indel_tbis.collect(),
        params.vqsr_indel_resource_labels
    )

    // Apply VQSR for INDELs
    GATK4_APPLYVQSR_INDEL(
        GATK4_APPLYVQSR_SNP.out.vcf
            .join(GATK4_APPLYVQSR_SNP.out.tbi)
            .combine(GATK4_VARIANTRECALIBRATOR_INDEL.out.recal.map { _meta, recal -> recal })
            .combine(GATK4_VARIANTRECALIBRATOR_INDEL.out.idx.map { _meta, idx -> idx })
            .combine(GATK4_VARIANTRECALIBRATOR_INDEL.out.tranches.map { _meta, tranches -> tranches }),
        ch_ref_fasta.map{ _meta, fasta -> fasta },
        ch_ref_fai.map{ _meta, fai -> fai },
        ch_ref_dict.map{ _meta, dict -> dict }
    )

    emit:
    vcf_filter_added = GATK4_APPLYVQSR_INDEL.out.vcf // [[meta], sorted.vcf.gz] - with filter tags
    vcf_filter_added_tbi = GATK4_APPLYVQSR_INDEL.out.tbi // [[meta], sorted.vcf.gz.tbi] - with filter tags
    // vcf_filtered         = // [[meta], sorted.vcf.gz] - filtered
    // vcf_filtered_tbi     = // [[meta], sorted.vcf.gz.tbi] - filtered
}
