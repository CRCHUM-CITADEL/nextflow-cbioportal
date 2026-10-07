process SIGPROFILER_DBS {
    tag "$meta.sample"
    label 'process_low'

    container params.container_sigprofiler

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(somatic_vcf)
        path signatures_db  // COSMIC v3.6 DBS-78 reference CSV
        path dbs_metadata

    output:
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_contribution_DBS.txt"), emit: sigs_dbs
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_counts_DBS.txt"), emit: sigs_counts_dbs

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    run_sigprofiler_dbs.py \\
        --vcf             ${somatic_vcf} \\
        --signatures_db   ${signatures_db} \\
        --metadata        ${dbs_metadata} \\
        --sample          ${meta.sample} \\
        --output_contrib  ${meta.sample}.data_mutational_signatures_contribution_DBS.txt \\
        --output_counts   ${meta.sample}.data_mutational_signatures_counts_DBS.txt
    """

    stub:
    """
    touch ${meta.sample}.data_mutational_signatures_contribution_DBS.txt
    touch ${meta.sample}.data_mutational_signatures_counts_DBS.txt
    """
}
