process FILTER_MUTATIONS {
    tag "${meta.sample_id}"
    label 'python'
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(maf), path(tsv)
    path(gene_panel)  // gene panel TSV (adds info.on_panel), or [] when unset

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    def exclude = params.exclude_genes && params.exclude_genes.toString() != 'null' ? params.exclude_genes.toString() : ''  // exclude_genes, '' = none
    """
    # Keep MAF rows overlapping a TSV row with Depth >= min_depth and VAF > min_vaf
    filter_mutations.py ${maf} "${tsv}" "${meta.sample_id}_mutations.txt" \\
        ${params.mutation_min_depth} ${params.mutation_min_vaf}
    filter_gene_panel.py "${meta.sample_id}_mutations.txt" --panel "${gene_panel}" --exclude "${exclude}" --type mutation --out "${meta.sample_id}_mutations.txt.tmp" && mv "${meta.sample_id}_mutations.txt.tmp" "${meta.sample_id}_mutations.txt"
    """
}
