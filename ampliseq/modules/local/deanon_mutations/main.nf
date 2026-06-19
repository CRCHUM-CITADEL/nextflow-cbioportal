process DEANON_MUTATIONS {
    label 'python'
    stageInMode 'copy'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(mutations_file)
    path(linking_file)

    output:
    path("data_mutations.txt")

    script:
    """
    format_mutations.py ${mutations_file} ${linking_file}
    """
}
