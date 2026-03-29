process MERGE_MUTATIONS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(mutation_files)

    output:
    path("data_mutations.txt")

    script:
    """
    files=( *_mutations.txt )
    head -1 "\${files[0]}" > data_mutations.txt
    for f in "\${files[@]}"; do
        tail -n +2 "\$f" >> data_mutations.txt
    done
    """
}
