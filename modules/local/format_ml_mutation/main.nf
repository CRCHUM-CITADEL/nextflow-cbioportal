process FORMAT_ML_MUTATION {
    publishDir "${params.outdir}/${group}/machine_learning/", mode: 'copy'
    
    container params.container_r

    input:
        val group
        path results_mutation

    output:
        path "all_somatic_mutations.tsv"
        path "all_germline_mutations.tsv"

    script:
    """
    ml_format_mutation.R $results_mutation
    """
}
