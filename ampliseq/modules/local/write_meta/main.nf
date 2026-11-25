process WRITE_META {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val(study_id)

    output:
    path("meta_*.txt")

    script:
    def no_seg = params.include_seg.toString().toBoolean() ? '' : '--no-seg'
    // Namespaces make cBioPortal show ANNOTATE_GNOMAD's gnomAD.* and ClinVar.* columns
    def mut_ns = params.filter_gnomad_af.toString().toBoolean() ? ['gnomAD', 'ClinVar'] : []
    // FusionMetrics.ZR_score only has values with a fusion-metrics TSV
    def sv_ns = params.fusion_metrics ? ['FusionMetrics'] : []
    def cna_ns = []
    // info.on_panel is only written with a gene panel (the validator rejects a namespace without columns)
    if (params.gene_panel && params.gene_panel.toString() != 'null') {
        mut_ns += 'info'; sv_ns += 'info'; cna_ns += 'info'
    }
    def namespaces = mut_ns ? "--namespaces ${mut_ns.join(',')}" : ''
    def sv_namespaces = sv_ns ? "--sv-namespaces ${sv_ns.join(',')}" : ''
    def cna_namespaces = cna_ns ? "--cna-namespaces ${cna_ns.join(',')}" : ''
    """
    format_meta.py "${study_id}" ${no_seg} ${namespaces} ${sv_namespaces} ${cna_namespaces}
    """
}
