process FORMAT_ML_CNV {
    publishDir "${params.outdir}/${group}/machine_learning/", mode: 'copy'
    
    container params.container_r

    input:
        val group
        path results_cnv_long

    output:
        path "cnv_gene.tsv"

    script:
    """
    ml_format_cnv.R $results_cnv_long
    """
}