process FORMAT_ML_MUTATION {
    publishDir "${params.outdir}/${group}/machine_learning/", mode: 'copy'

    container params.container_r

    input:
<<<<<<< HEAD
        tuple val(group), path(results_mutation)
=======
        val group
        path results_mutation
>>>>>>> 555063098b4c3191c53c1bb5f99da149294529b1

    output:
        path "all_somatic_mutations.tsv"
        path "all_germline_mutations.tsv"

    script:
    """
    ml_format_mutation.R $results_mutation
    """
}
