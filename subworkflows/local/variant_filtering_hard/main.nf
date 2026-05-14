//
// SUBWORKFLOW: VARIANT_FILTERING_HARD
//
// Hard-filter germline short variants (SNPs + indels separately, then recombine)
// Utilises scatter-gather approach, i.e. input should be per interval
//

include { GATK4_SELECTVARIANTS as SELECT_SNP              } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as SELECT_INDEL            } from '../../../modules/nf-core/gatk4/selectvariants/main'
include { GATK4_SELECTVARIANTS as EXCLUDE_FILTERED_SNP    } from '../../../modules/nf-core/gatk4/selectvariants/'
include { GATK4_SELECTVARIANTS as EXCLUDE_FILTERED_INDEL  } from '../../../modules/nf-core/gatk4/selectvariants/'
include { GATK4_VARIANTFILTRATION as FILTER_SNP           } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_VARIANTFILTRATION as FILTER_INDEL         } from '../../../modules/nf-core/gatk4/variantfiltration/main'
include { GATK4_MERGEVCFS as MERGEVCFS_FILTER_ADDED       } from '../../../modules/nf-core/gatk4/mergevcfs/main'
include { GATK4_MERGEVCFS as MERGEVCFS_FILTERED           } from '../../../modules/nf-core/gatk4/mergevcfs/main'

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
    MERGEVCFS_FILTER_ADDED(
        ch_filter_added,   // [[meta], [snp.vcf.gz, indel.vcf.gz]]
        ch_ref_dict,
    )

    // Rejoin SNP and indels with filters removed
    def ch_filtered = EXCLUDE_FILTERED_SNP.out.vcf
        .join(EXCLUDE_FILTERED_INDEL.out.vcf)
        .map { meta, snp_vcf, indel_vcf -> [meta, [snp_vcf, indel_vcf]] }
    MERGEVCFS_FILTERED(
        ch_filtered,       // [[meta], [snp.vcf.gz, indel.vcf.gz]]
        ch_ref_dict,
    )

    emit:
    vcf_filter_added     = MERGEVCFS_FILTER_ADDED.out.vcf   // [[meta], sorted.vcf.gz] - with filter tags
    vcf_filter_added_tbi = MERGEVCFS_FILTER_ADDED.out.tbi   // [[meta], sorted.vcf.gz.tbi] - filter tags
    vcf_filtered         = MERGEVCFS_FILTERED.out.vcf       // [[meta], sorted.vcf.gz] - filtered
    vcf_filtered_tbi     = MERGEVCFS_FILTERED.out.tbi       // [[meta], sorted.vcf.gz.tbi] - filtered
}
