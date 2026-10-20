process FILTER_MUTATIONS {
    tag "${meta.sample_id}"
    label 'python'
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(maf), path(tsv)

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    """
    # Keep MAF rows overlapping a TSV row with Depth >= min_depth and VAF > min_vaf
    filter_mutations.py ${maf} "${tsv}" "${meta.sample_id}_mutations.txt" \\
        ${params.mutation_min_depth} ${params.mutation_min_vaf}
    """
}
