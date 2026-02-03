process FORMAT_ML_CNV {
    publishDir "${params.outdir}/${group}/machine_learning/", mode: 'copy'

    container params.container_r

    input:
        tuple val(group), path(mutation_results)

    output:
        path "mutations_processes_*.tsv"

    script:
    """
    ml_format_cnv.R $results_cnv_long 
    """
}
