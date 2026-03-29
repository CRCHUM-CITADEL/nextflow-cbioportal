process WRITE_META {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val(study_id)

    output:
    path("meta_*.txt")

    script:
    """
    format_meta.py "${study_id}"
    """
}
