process FORMAT_ML_EXPRESSION {
    publishDir "${params.outdir}/${group}/machine_learning/formatted", mode: 'copy'

    container params.container_r

    input:
        tuple val(group), path(results_expression)

    output:
        tuple val(group), path("expression_tpm.tsv")

    script:
    """
    ml_format_expression.R $results_expression
    """

    stub:
    """
    touch expression_tpm.tsv
    """
}
