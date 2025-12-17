process FORMAT_ML_EXPRESSION {
    publishDir "${params.outdir}/${group}/ml/", mode: 'copy'

    container params.container_r

    input:
        val group
        path results_expression

    output:
        path "expression_tpm.tsv"

    script:
    """
    ml_format_expression.R $results_expression
    """
}
