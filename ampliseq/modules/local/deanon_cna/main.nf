process DEANON_CNA {
    label 'python'
    stageInMode 'copy'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(cna_file)
    path(linking_file)

    output:
    path("data_cna.txt")

    script:
    """
    format_cna_deanon.py ${cna_file} ${linking_file}
    """
}
