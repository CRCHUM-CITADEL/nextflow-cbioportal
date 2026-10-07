process SIGPROFILER_ID {
    tag "$meta.sample"
    label 'process_low'

    container params.container_sigprofiler

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(somatic_vcf)
        path signatures_db  // COSMIC v3.6 ID-83 reference CSV
        path id_metadata
        path fasta

    output:
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_contribution_ID.txt"), emit: sigs_id
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_counts_ID.txt"), emit: sigs_counts_id

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    run_sigprofiler_id.py \\
        --vcf             ${somatic_vcf} \\
        --fasta           ${fasta} \\
        --signatures_db   ${signatures_db} \\
        --metadata        ${id_metadata} \\
        --sample          ${meta.sample} \\
        --output_contrib  ${meta.sample}.data_mutational_signatures_contribution_ID.txt \\
        --output_counts   ${meta.sample}.data_mutational_signatures_counts_ID.txt
    """

    stub:
    """
    touch ${meta.sample}.data_mutational_signatures_contribution_ID.txt
    touch ${meta.sample}.data_mutational_signatures_counts_ID.txt
    """
}
