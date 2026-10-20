process WRITE_META {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val(study_id)

    output:
    path("meta_*.txt")

    script:
    def no_seg = params.include_seg.toString().toBoolean() ? '' : '--no-seg'
    // gnomAD namespace: cBioPortal shows ANNOTATE_GNOMAD's gnomAD.AF column
    def namespaces = params.filter_gnomad_af.toString().toBoolean() ? '--namespaces gnomAD' : ''
    """
    format_meta.py "${study_id}" ${no_seg} ${namespaces}
    """
}
