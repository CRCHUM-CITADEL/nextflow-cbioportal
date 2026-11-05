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
    def namespaces = params.filter_gnomad_af.toString().toBoolean() ? '--namespaces gnomAD,ClinVar' : ''
    // FusionMetrics.ZR_score only has values with a fusion-metrics TSV
    def sv_namespaces = params.fusion_metrics ? '--sv-namespaces FusionMetrics' : ''
    """
    format_meta.py "${study_id}" ${no_seg} ${namespaces} ${sv_namespaces}
    """
}
