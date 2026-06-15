process GATK4_VARIANTRECALIBRATOR {
    // This process builds a VQSR model using GATK VariantRecalibrator, based on
    // input VCFs and parameters, and outputs a tarball containing the recalibration
    // model, tranches, R plots, and related files for later variant quality score
    // recalibration.
    //
    // Based on MalariaGEN Pf8 implementation - https://github.com/malariagen/malariagen-pf8-snp-indel-calling/blob/master/modules/variant_recalibration_build_model.nf (MIT license)
    // and the existing nf-core module - https://nf-co.re/modules/gatk4_variantrecalibrator/
    tag "${meta.id}_${mode}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/ce/ced519873646379e287bc28738bdf88e975edd39a92e7bc6a34bccd37153d9d0/data'
        : 'community.wave.seqera.io/library/gatk4_gcnvkernel:edb12e4f0bf02cd3'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    val mode
    path fasta
    path fai
    path dict
    path resource_vcfs
    path resource_tbis
    val resource_labels

    output:
    tuple val(meta), path("*.recal"),       emit: recal
    tuple val(meta), path("*.recal.idx"),   emit: idx
    tuple val(meta), path("*.tranches"),    emit: tranches
    tuple val(meta), path("*plots.R"), emit: plots, optional: true
    tuple val("${task.process}"), val('gatk4'), eval("gatk --version | sed -n '/GATK.*v/s/.*v//p'"), topic: versions, emit: versions_gatk4

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${mode.toLowerCase()}"

    // Memory settings adapter from nf-core GATK4_VARIANTRECALIBRATOR
    // Alternatively, the MalariaGEN Pf8 approach could be used: https://github.com/malariagen/malariagen-pf8-snp-indel-calling/blob/529fe6b59bbf14a99d8c46264f6a38c4761a0ffa/modules/variant_recalibration_build_model.nf#L20
    def avail_mem = 3072
    if (!task.memory) {
        log.info('[GATK VariantRecalibrator] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.')
    } else {
        avail_mem = (task.memory.mega * 0.8).intValue()
    }

    // Sanity checking provided resource inputs and labels
    if (resource_labels.size() != resource_vcfs.size()) {
        error "Number of resource labels (${resource_labels.size()}) must match number of resource VCFs (${resource_vcfs.size()})"
    }

    // Build resource commands
    def resource_commands = []
    resource_labels.eachWithIndex { label, index ->
        if (index < resource_vcfs.size()) {
            resource_commands << "${label} ${resource_vcfs[index]}"
        }
    }

    """
    gatk --java-options "-Xmx${avail_mem}M -XX:-UsePerfData" \\
        VariantRecalibrator \\
        --variant ${vcf} \\
        --output ${prefix}.recal \\
        --tranches-file ${prefix}.tranches \\
        --mode ${mode} \\
        --reference ${fasta} \\
        ${resource_commands.join(' ')} \\
        --tmp-dir . \\
        ${args}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}_${mode.toLowerCase()}"
    """
    touch ${prefix}.recal
    touch ${prefix}.recal.idx
    touch ${prefix}.tranches
    touch ${prefix}plots.R
    """
}
