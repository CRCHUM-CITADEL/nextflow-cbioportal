process PROCESS_ML_EXPRESSION {
    publishDir "${params.outdir}/${group}/machine_learning/processed/", mode: 'copy'

    container params.container_r

    input:
        tuple val(group), path(results_expression)

    output:
        path "expression_processed_log2.tsv", emit: log2
        path "expression_processed_standardized.tsv", emit: std

    script:
    """
    ml_expression_processor.R $results_expression
    """
}
