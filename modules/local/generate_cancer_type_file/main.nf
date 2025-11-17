process GENERATE_CANCER_TYPE_FILE { 
    publishDir "${params.outdir}/${group}", mode: 'copy'

    input:
    val group

    output:
    path "cancer_type.txt"

    script:
    """
    echo -e "${group}\tPlaceholder cancer type\tOrangeRed\ttissue" > cancer_type.txt
    """
}
