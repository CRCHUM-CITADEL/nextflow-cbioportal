process PROCESS_ML_MUTATION {
    publishDir "${params.outdir}/${group}/machine_learning/processed", mode: 'copy'

    container params.container_r

    input:
        tuple val(group), path(mutation_results)

    output:
        path "mutations_processed_*.tsv"

    script:
    """
    ml_mutation_processor.R $mutation_results
    """
}
