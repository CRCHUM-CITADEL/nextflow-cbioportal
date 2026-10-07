process PROCESS_ML_CNV {
    publishDir { "${params.outdir}/${group}/machine_learning/processed/" }, mode: 'copy'

    label "process_high_memory"
    container params.container_r

    input:
        tuple val(group), path(results_cnv)

    output:
        path "cnv_genes.tsv", emit: tsv

    script:
    """
    ml_cnv_processor.R $results_cnv
    """

    stub:
    """
    touch cnv_genes.tsv
    """
}
