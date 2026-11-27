process FORMAT_CNA {
    tag "${meta.sample_id}"
    label 'python'
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(tsv)
    path(gene_panel)  // gene panel TSV (adds info.on_panel), or [] when unset

    output:
    path("${meta.sample_id}_cna.txt")

    script:
    def exclude = params.exclude_genes && params.exclude_genes.toString() != 'null' ? params.exclude_genes.toString() : ''  // exclude_genes, '' = none
    // Filter mode: CN >= cna_min_copy_number only (toBoolean: NF 26 CLI params are strings)
    def filter_args = params.filter_tsv_variants.toString().toBoolean()
        ? "--min-copy-number ${params.cna_min_copy_number}"
        : ''
    """
    format_cna.py "${tsv}" "${meta.sample_id}" ${filter_args}
    mv data_cna.txt "${meta.sample_id}_cna.txt"
    filter_gene_panel.py "${meta.sample_id}_cna.txt" --panel "${gene_panel}" --exclude "${exclude}" --type cna --out "${meta.sample_id}_cna.txt.tmp" && mv "${meta.sample_id}_cna.txt.tmp" "${meta.sample_id}_cna.txt"
    """
}
