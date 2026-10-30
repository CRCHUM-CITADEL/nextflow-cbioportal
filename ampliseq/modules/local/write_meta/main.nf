process WRITE_META {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val(study_id)

    output:
    path("meta_*.txt")

    script:
    def no_seg = params.include_seg.toString().toBoolean() ? '' : '--no-seg'
    // Namespaces make cBioPortal show ANNOTATE_GNOMAD's gnomAD.* and ClinVar.* columns; either filter runs it
    def annotated  = params.filter_gnomad_af.toString().toBoolean() || params.filter_clinvar.toString().toBoolean()
    def namespaces = annotated ? '--namespaces gnomAD,ClinVar' : ''
    """
    format_meta.py "${study_id}" ${no_seg} ${namespaces}
    """
}
