process SNPEFF_BUILD {
    tag "$db_name"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/30/30669e5208952f30d59d0d559928772f082830d01a140a853fff13a2283a17b0/data'
        : 'community.wave.seqera.io/library/snpeff:5.4.0a--eaf6ce30125b2b17'}"

    input:
    tuple val(meta_ref), path(fasta), path(annotation), path(cds)
    val db_name
    val annotation_format           // 'gff', 'gtf' or empty (falls back to detection in filename)
    tuple val(meta_bed), path(bed)  // BED file for mitochondrial/apicoplast detection
    path snpeff_config_template

    output:
    tuple val(meta_ref), path("snpeff_db"), emit: db
    tuple val("${task.process}"), val('snpeff'), eval("snpEff -version 2>&1 | cut -f 2 -d '\t'"), topic: versions, emit: versions_snpeff

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    // extract gff/gtf format from filename if not provided
    if (!annotation_format) {
        def anno_name = annotation.name.toLowerCase()
        if (anno_name.endsWith('.gtf') || anno_name.endsWith('.gtf.gz')) {
            annotation_format = 'gtf'
        } else if (anno_name.endsWith('.gff') || anno_name.endsWith('.gff.gz') ||
                anno_name.endsWith('.gff3') || anno_name.endsWith('.gff3.gz')) {
            annotation_format = 'gff'
        } else {
            error "Could not determine annotation format from filename: ${annotation.name}. " +
                "Please provide annotation_format parameter ('gtf' or 'gff')."
        }
    }
    if (annotation_format != 'gtf' && annotation_format != 'gff') {
        error "Invalid annotation_format: '${annotation_format}'. Must be 'gtf' or 'gff'."
    }
    def annotation_file = (annotation_format == 'gtf') ? 'genes.gtf' : 'genes.gff'
    def annotation_arg = (annotation_format == 'gtf') ? '-gtf22' : '-gff3'
    """
    # Create the directory structure snpEff expects
    mkdir -p snpeff_db/data/${db_name}

    # Copy and rename files to match snpEff naming requirements, unzipping them if necessary
    if [[ ${fasta} == *.gz ]]; then
        gunzip -c ${fasta} > snpeff_db/data/${db_name}/sequences.fa
    else
        cp ${fasta} snpeff_db/data/${db_name}/sequences.fa
    fi
    if [[ ${annotation} == *.gz ]]; then
        gunzip -c ${annotation} > snpeff_db/data/${db_name}/${annotation_file}
    else
        cp ${annotation} snpeff_db/data/${db_name}/${annotation_file}
    fi
    if [[ ${cds} == *.gz ]]; then
        gunzip -c ${cds} > snpeff_db/data/${db_name}/cds.fa
    else
        cp ${cds} snpeff_db/data/${db_name}/cds.fa
    fi

    # Create snpEff config file, starting with the template file
    cp ${snpeff_config_template} snpeff_db/snpEff.config

    # Append custom genome configuration to config
    cat >> snpeff_db/snpEff.config << EOF

# ${db_name} genome configuration
${db_name}.genome : ${db_name}
EOF

    # Extract mitochondrial and apicoplast contig names from BED file (case-insensitive)
    MITO_CONTIGS=\$(awk '\$1 ~ /[Mm][Ii][Tt]/ {print \$1}' ${bed} | sort -u)
    API_CONTIGS=\$(awk '\$1 ~ /[Aa][Pp][Ii]/ {print \$1}' ${bed} | sort -u)

    # Add mito and api codon tables if found
    if [ -n "\${MITO_CONTIGS}" ]; then
        echo "# Mitochondrial codon table configuration" >> snpeff_db/snpEff.config
        for contig in \${MITO_CONTIGS}; do
            echo "    ${db_name}.\${contig}.codonTable : Protozoan_Mitochondrial" >> snpeff_db/snpEff.config
        done
    fi
    if [ -n "\${API_CONTIGS}" ]; then
        echo "# Apicoplast codon table configuration" >> snpeff_db/snpEff.config
        for contig in \${API_CONTIGS}; do
            echo "    ${db_name}.\${contig}.codonTable : Bacterial_and_Plant_Plastid" >> snpeff_db/snpEff.config
        done
    fi

    # Build the database - dataDir is relative to the location of snpEff.config
    snpEff build \\
        -c snpeff_db/snpEff.config \\
        ${annotation_arg} \\
        -v \\
        -noCheckProtein \\
        -dataDir ./data/ \\
        ${args} \\
        ${db_name}
    """

    stub:
    // Handle empty annotation_format in stub
    def stub_format = annotation_format ?: 'gtf'
    def annotation_file = (stub_format == 'gtf') ? 'genes.gtf' : 'genes.gff'
    """
    mkdir -p snpeff_db/data/${db_name}
    touch snpeff_db/data/${db_name}/sequences.fa
    touch snpeff_db/data/${db_name}/${annotation_file}
    touch snpeff_db/data/${db_name}/cds.fa
    touch snpeff_db/snpEff.config
    touch versions.yml
    """
}
