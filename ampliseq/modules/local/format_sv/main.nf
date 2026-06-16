process FORMAT_SV {
    tag "${meta.sample_id}"
    label 'python'
    publishDir "${params.outdir}/samples/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(tsv)

    output:
    path("${meta.sample_id}_sv.txt")

    script:
    """
    format_tsv.py "${tsv}" "${meta.subject_id}"
    mv data_sv.txt "${meta.sample_id}_sv.txt"
    """
}
