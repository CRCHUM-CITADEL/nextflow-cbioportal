process FORMAT_ML_MUTATION {
    publishDir { "${params.outdir}/${group}/machine_learning/formatted" }, mode: 'copy'

    container params.container_r
    label "process_medium_memory"
    input:
        tuple val(group), path(results_mutation)

    output:
        tuple val(group), path("all_somatic_mutations.tsv")
        // TODO : est-ce qu'on mets germline?
        // tuple val(group), path("all_germline_mutations.tsv")

    script:
    """
    ml_format_mutation.R $results_mutation
    """

    stub:
    """
    touch all_somatic_mutations.tsv
    """
}
