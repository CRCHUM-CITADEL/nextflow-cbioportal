process PASSTHROUGH_MUTATIONS {
    tag "${meta.sample_id}"
    label 'python'
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(maf)
    path(gene_panel)  // gene panel TSV (adds info.on_panel), or [] when unset

    output:
    path("${meta.sample_id}_mutations.txt")

    script:
    def exclude = params.exclude_genes && params.exclude_genes.toString() != 'null' ? params.exclude_genes.toString() : ''  // exclude_genes, '' = none
    """
    cp "${maf}" "${meta.sample_id}_mutations.txt"
    filter_gene_panel.py "${meta.sample_id}_mutations.txt" --panel "${gene_panel}" --exclude "${exclude}" --type mutation --out "${meta.sample_id}_mutations.txt.tmp" && mv "${meta.sample_id}_mutations.txt.tmp" "${meta.sample_id}_mutations.txt"
    """
}
