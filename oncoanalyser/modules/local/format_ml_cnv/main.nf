process FORMAT_ML_CNV {
    publishDir "${params.outdir}/${group}/machine_learning/formatted", mode: 'copy'

    label "process_high_memory"

    container params.container_r

    input:
        tuple val(group), path(results_cnv_long)

    output:
        tuple val(group), path("cnv_gene.tsv")

    script:
    """
    ml_format_cnv.R $results_cnv_long
    """

    stub:
    """
    touch cnv_gene.tsv
    """
}
