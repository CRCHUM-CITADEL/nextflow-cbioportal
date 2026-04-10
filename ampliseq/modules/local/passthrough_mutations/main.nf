process PASSTHROUGH_MUTATIONS {
    tag "${meta.sample_id}"
    publishDir "${params.outdir}/samples/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(maf)

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    """
    cp "${maf}" "${meta.sample_id}_mutations.txt"
    """
}
