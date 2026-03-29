process FORMAT_SV {
    tag "${meta.sample_id}"
    label 'python'

    input:
    tuple val(meta), path(tsv)

    output:
    path("${meta.sample_id}_sv.txt")

    script:
    """
    format_tsv.py "${tsv}" "${meta.sample_id}"
    mv data_sv.txt "${meta.sample_id}_sv.txt"
    """
}
