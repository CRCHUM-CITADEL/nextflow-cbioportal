process WRITE_GENE_PANELS {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(gene_panel)
    path(linking_file)
    val(cnv_failed_ids) // sample IDs that failed CNV depth QC: cna = NA in the matrix (may be empty)
    val(study_id)

    output:
    path("data_gene_panel_matrix.txt"), emit: matrix
    path("meta_gene_panel_matrix.txt"), emit: meta
    path("gene_panels"),                emit: panels

    script:
    def failed = cnv_failed_ids ? "--cnv-failed ${cnv_failed_ids.join(' ')}" : ''
    """
    write_gene_panels.py ${gene_panel} ${linking_file} "${study_id}" ${failed}
    """

    stub:
    """
    mkdir -p gene_panels
    touch data_gene_panel_matrix.txt meta_gene_panel_matrix.txt gene_panels/gene_panel_${study_id}_mutation.txt
    """
}
