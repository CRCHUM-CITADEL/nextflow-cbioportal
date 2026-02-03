process FORMAT_ML_EXPRESSION {
    publishDir "${params.outdir}/${group}/machine_learning/", mode: 'copy'

    container params.container_r

    input:
<<<<<<< HEAD
        tuple val(group), path(results_expression)
=======
        val group
        path results_expression
>>>>>>> 555063098b4c3191c53c1bb5f99da149294529b1

    output:
        path "expression_tpm.tsv"

    script:
    """
    ml_format_expression.R $results_expression
    """
}
