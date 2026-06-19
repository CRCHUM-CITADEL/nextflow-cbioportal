process MERGE_SV {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(sv_files)

    output:
    path("data_sv.txt")

    script:
    """
    files=( *_sv.txt )
    head -1 "\${files[0]}" > data_sv.txt
    for f in "\${files[@]}"; do
        tail -n +2 "\$f" >> data_sv.txt
    done
    """
}
