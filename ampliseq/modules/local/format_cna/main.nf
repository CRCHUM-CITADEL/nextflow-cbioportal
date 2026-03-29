process FORMAT_CNA {
    tag "${meta.sample_id}"
    label 'python'

    input:
    tuple val(meta), path(tsv)

    output:
    path("${meta.sample_id}_cna.txt")

    script:
    """
    format_cna.py "${tsv}" "${meta.sample_id}"
    mv data_cna.txt "${meta.sample_id}_cna.txt"
    """
}
