process GENERATE_META_FILE {
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    val label
    val text

    output:
    path "meta_${label}.txt"

    script:
    """
    echo -e "${text}" > meta_${label}.txt
    """
}
