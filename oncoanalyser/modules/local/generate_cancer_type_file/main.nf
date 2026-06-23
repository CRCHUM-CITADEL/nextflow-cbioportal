process GENERATE_CANCER_TYPE_FILE {
    publishDir "${params.outdir}/${group}", mode: 'copy'

    input:
    val group
    val oncotree_code

    output:
    path "cancer_type.txt"

    script:
    if (!oncotree_code) error "params.oncotree_code is required for clinical mode. Set it to a valid OncoTree code (e.g., 'CHOL'). See https://oncotree.info"
    def type_id = oncotree_code.toLowerCase()
    """
    echo -e "${type_id}\t${oncotree_code}\tOrangeRed\ttissue" > cancer_type.txt
    """
}
