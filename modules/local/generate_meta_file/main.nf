process GENERATE_META_FILE {
    publishDir "${params.outdir}/${group}", mode: 'copy'

    input:
    val group
    val label
    val text

    output:
    path "meta_${label}.txt"

    script:
    """
sed "s/add_text/${group}/g" << EOF > identified_text.txt
${text}
EOF
    cat identified_text.txt > meta_${label}.txt
    """
}
