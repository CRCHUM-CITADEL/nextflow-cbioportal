process MERGE_CNA {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(cna_files)

    output:
    path("data_cna.txt")

    script:
    """
    files=( *_cna.txt )
    head -1 "\${files[0]}" > data_cna.txt
    for f in "\${files[@]}"; do
        tail -n +2 "\$f" >> data_cna.txt
    done
    """
}
