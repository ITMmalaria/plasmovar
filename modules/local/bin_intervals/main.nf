process BIN_INTERVALS {
    tag "${meta.id}"
    label 'process_single'

    conda "conda-forge::python=3.12"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.12' :
        'biocontainers/python:3.12' }"

    input:
    tuple val(meta), path(bed)
    val min_contig_size
    val target_interval_size

    output:
    path "*.bed", emit: intervals
    path "*.txt", emit: stats

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def target_bin_size = target_interval_size ? "--target-bin-size ${target_interval_size}" : ""
    """
    prepare_intervals.py \\
        --bed ${bed} \\
        --min-contig-size ${min_contig_size} \\
        ${target_bin_size} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
