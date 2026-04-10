process DEANON_SEG {
    label 'python'
    stageInMode 'copy'
    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
    path(seg_file)
    path(linking_file)

    output:
    path("data_seg.txt")

    script:
    """
    seg_deanon.py ${seg_file} ${linking_file}
    """
}
