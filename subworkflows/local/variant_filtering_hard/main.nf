//
// SUBWORKFLOW: VARIANT_FILTERING_HARD
//
// Hard-filter germline short variants (SNPs + indels separately, then recombine)
// Utilises scatter-gather approach, i.e. input should be per interval
//

include { GATK4_SELECTVARIANTS as SELECT_SNP                     } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_INDEL                   } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as EXCLUDE_FILTERED_SNP           } from '../../../modules/nf-core/gatk4/selectvariants/'
include { GATK4_SELECTVARIANTS as EXCLUDE_FILTERED_INDEL         } from '../../../modules/nf-core/gatk4/selectvariants/'
include { GATK4_VARIANTFILTRATION as FILTER_SNP                  } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION as FILTER_INDEL                } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_MERGEVCFS as CONCAT_VCFS_SNP_INDEL_FILTER_ADDED   } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as CONCAT_VCFS_SNP_INDEL_FILTERED       } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as GATHER_VCFS_BY_INTERVAL_FILTER_ADDED } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as GATHER_VCFS_BY_INTERVAL_FILTERED     } from '../../../modules/nf-core/gatk4/mergevcfs/main'


workflow VARIANT_FILTERING_HARD {

    take:
    ch_vcf_by_interval  // [[meta], vcf.gz, vcf.gz.tbi] - 1 per interval
    ch_ref_fasta        // [[meta], fasta]
    ch_ref_fai          // [[meta], fasta.fai]
    ch_ref_dict         // [[meta], dict]

    main:

    // Fork SNPs and indels so they can be processed separately
    SELECT_SNP(
        ch_vcf_by_interval.map { meta, vcf, tbi -> [meta, vcf, tbi, []] },
    )
    SELECT_INDEL(
        ch_vcf_by_interval.map { meta, vcf, tbi -> [meta, vcf, tbi, []] },
    )

    // Annotate VCFs with filter tags
    FILTER_SNP(
        SELECT_SNP.out.vcf.join(SELECT_SNP.out.tbi),
        ch_ref_fasta,
        ch_ref_fai,
        ch_ref_dict,
        [[:], []],  // gzi index is only needed if fasta input is in bgzip format
        // ch_gzi,
    )
    FILTER_INDEL(
        SELECT_INDEL.out.vcf.join(SELECT_INDEL.out.tbi),
        ch_ref_fasta,
        ch_ref_fai,
        ch_ref_dict,
        [[:], []],  // ch_gzi,
    )

    // Exclude filtered variants from VCF
    EXCLUDE_FILTERED_SNP(
        FILTER_SNP.out.vcf
            .join(FILTER_SNP.out.tbi)
            .map { meta, vcf, tbi -> [meta, vcf, tbi, []] },
    )
    EXCLUDE_FILTERED_INDEL(
        FILTER_INDEL.out.vcf
            .join(FILTER_INDEL.out.tbi)
            .map { meta, vcf, tbi -> [meta, vcf, tbi, []] },
    )

    // Rejoin SNP and indels annotated with FILTER column
    def ch_filter_added = FILTER_SNP.out.vcf
        .join(FILTER_INDEL.out.vcf)
        .map { meta, snp_vcf, indel_vcf -> [meta, [snp_vcf, indel_vcf]] }
    CONCAT_VCFS_SNP_INDEL_FILTER_ADDED(
        ch_filter_added,   // [[meta], [snp.vcf.gz, indel.vcf.gz]]
        ch_ref_dict,
    )

    // Rejoin SNP and indels with filters removed
    def ch_filtered = EXCLUDE_FILTERED_SNP.out.vcf
        .join(EXCLUDE_FILTERED_INDEL.out.vcf)
        .map { meta, snp_vcf, indel_vcf -> [meta, [snp_vcf, indel_vcf]] }
    CONCAT_VCFS_SNP_INDEL_FILTERED(
        ch_filtered,       // [[meta], [snp.vcf.gz, indel.vcf.gz]]
        ch_ref_dict,
    )

    // Prepare for VCF merging: gather all interval VCFs
    // Group by genome_id to merge all intervals for that genome
    ch_vcfs_to_merge_filter_added = CONCAT_VCFS_SNP_INDEL_FILTER_ADDED.out.vcf
        // [ meta, vcf ]    (1 per interval)
        .map { meta, vcf ->
            [ meta.genome_id, vcf ]
        }
        .groupTuple(by: 0)
        // [ meta.genome_id, [vcfs] ] (1 element with all interval files)
        .map { genome_id, vcfs ->
            def final_meta = [ id: "${genome_id}_merged", genome_id: genome_id ]
            [ final_meta, vcfs ]
        }
    ch_vcfs_to_merge_filtered = CONCAT_VCFS_SNP_INDEL_FILTERED.out.vcf
        // [ meta, vcf ]    (1 per interval)
        .map { meta, vcf ->
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
    GATHER_VCFS_BY_INTERVAL_FILTER_ADDED(
        ch_vcfs_to_merge_filter_added,
        ch_ref_dict
    )
    GATHER_VCFS_BY_INTERVAL_FILTERED(
        ch_vcfs_to_merge_filtered,
        ch_ref_dict
    )

    emit:
    vcf_filter_added     = GATHER_VCFS_BY_INTERVAL_FILTER_ADDED.out.vcf   // [[meta], sorted.vcf.gz] - with filter tags
    vcf_filter_added_tbi = GATHER_VCFS_BY_INTERVAL_FILTER_ADDED.out.tbi   // [[meta], sorted.vcf.gz.tbi] - filter tags
    vcf_filtered         = GATHER_VCFS_BY_INTERVAL_FILTERED.out.vcf       // [[meta], sorted.vcf.gz] - filtered
    vcf_filtered_tbi     = GATHER_VCFS_BY_INTERVAL_FILTERED.out.tbi       // [[meta], sorted.vcf.gz.tbi] - filtered
}
