process FILTER_GENE_PANEL {
    tag "${data_type}"
    label 'python'

    input:
    path(data_file)
    path(gene_panel)
    val(data_type)    // mutation, cna or sv

    output:
    path("filtered/${data_file.name}")

    script:
    """
    mkdir -p filtered
    filter_gene_panel.py ${data_file} --panel ${gene_panel} --type ${data_type} --out "filtered/${data_file.name}"
    """

    stub:
    """
    mkdir -p filtered
    cp ${data_file} "filtered/${data_file.name}"
    """
}
