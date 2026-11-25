process FILTER_GENE_PANEL {
    tag "${data_type}"
    label 'python'

    input:
    path(data_file)
    path(gene_panel)  // gene panel TSV (adds info.on_panel), or [] for exclude_genes only
    val(data_type)    // mutation, cna or sv

    output:
    path("filtered/${data_file.name}")

    script:
    def exclude = params.exclude_genes && params.exclude_genes.toString() != 'null' ? params.exclude_genes.toString() : ''  // exclude_genes, '' = none
    """
    mkdir -p filtered
    filter_gene_panel.py ${data_file} --panel "${gene_panel}" --exclude "${exclude}" --type ${data_type} --out "filtered/${data_file.name}"
    """

    stub:
    """
    mkdir -p filtered
    cp ${data_file} "filtered/${data_file.name}"
    """
}
