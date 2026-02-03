process FORMAT_ML_CNV {
    publishDir "${params.outdir}/${group}/machine_learning/", mode: 'copy'

    container params.container_r

    input:
<<<<<<< HEAD
        tuple val(group), path(results_cnv_long)
=======
        val group
        path results_cnv_long
>>>>>>> 555063098b4c3191c53c1bb5f99da149294529b1

    output:
        path "cnv_gene.tsv"

    script:
    """
    ml_format_cnv.R $results_cnv_long
    """
}
