process PASSTHROUGH_MUTATIONS {
    tag "${meta.sample_id}"

    input:
    tuple val(meta), path(maf)

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    """
    cp "${maf}" "${meta.sample_id}_mutations.txt"
    """
}
