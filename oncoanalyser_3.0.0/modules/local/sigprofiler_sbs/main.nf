process SIGPROFILER_SBS {
    tag "$meta.sample"
    label 'process_low'

    container params.container_sigprofiler

    publishDir { "${params.outdir}/${meta.group}/${meta.subject}" }, mode: 'copy'

    input:
        tuple val(meta), path(snv_counts)
        path signatures_db  // COSMIC v3.6 SBS-96 reference CSV
        path sbs_metadata

    output:
        tuple val(meta), path("${meta.sample}.data_mutational_signatures_contribution_SBS.txt"), emit: sigs

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    run_sigprofiler_sbs.py \\
        --snv_counts    ${snv_counts} \\
        --signatures_db ${signatures_db} \\
        --metadata      ${sbs_metadata} \\
        --sample        ${meta.sample} \\
        --output        ${meta.sample}.data_mutational_signatures_contribution_SBS.txt
    """

    stub:
    """
    touch ${meta.sample}.data_mutational_signatures_contribution_SBS.txt
    """
}
