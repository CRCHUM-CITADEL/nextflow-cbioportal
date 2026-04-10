process DEANON_SV {
    label 'python'
    stageInMode 'copy'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(sv_file)
    path(linking_file)

    output:
    path("data_sv.txt")

    script:
    """
    format_sv.py ${sv_file} ${linking_file}
    """
}
