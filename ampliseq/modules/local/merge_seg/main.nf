process MERGE_SEG {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(seg_files)

    output:
    path("data_seg.txt")

    script:
    """
    files=( *_seg.txt )
    head -1 "\${files[0]}" > data_seg.txt
    for f in "\${files[@]}"; do
        tail -n +2 "\$f" >> data_seg.txt
    done
    """
}
