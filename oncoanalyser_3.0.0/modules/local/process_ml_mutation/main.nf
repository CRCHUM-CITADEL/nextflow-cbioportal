process PROCESS_ML_MUTATION {
    publishDir { "${params.outdir}/${group}/machine_learning/processed" }, mode: 'copy'
    label "process_medium"

    container params.container_r

    input:
        tuple val(group), path(mutation_results)
        path hotspots_json  // pre-staged cancerhotspots.org data, or [] to fetch live

    output:
        path "mutations_processed_*.tsv"

    script:
    // hotspots_json is [] (falsy in Groovy) when params.cancer_hotspots_data is unset;
    // the R script falls back to a live GET when it gets no second argument.
    def hotspots_arg = hotspots_json ? "${hotspots_json}" : ""
    """
    ml_mutation_processor.R $mutation_results $hotspots_arg
    """

    stub:
    """
    touch mutations_processed_somatic.tsv
    """
}
