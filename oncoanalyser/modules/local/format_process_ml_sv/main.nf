process FORMAT_PROCESS_ML_SV {
    publishDir "${params.outdir}/${group}/machine_learning/processed", mode: 'copy'

    container params.container_r

    input:
        path cosmic_data
        path known_fusions
        tuple val(group), path(results_sv)

    output:
        path "filtered_fusions_matrix.tsv"

    script:
    """
    ml_format_sv_processor.R $results_sv $known_fusions $cosmic_data
    """
}
