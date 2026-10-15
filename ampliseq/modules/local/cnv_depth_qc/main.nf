process CNV_DEPTH_QC {
    tag "${meta.sample_id}"
    label 'python'

    input:
    tuple val(meta), path(sample_folder)

    output:
    path("${meta.sample_id}_cnv_depth_qc.txt")

    script:
    def pass_only = params.cnv_depth_pass_only.toString().toBoolean() ? '--pass-only' : ''
    """
    VCF=\$(find -L "${sample_folder}" -maxdepth 1 \\( -name '*-basespace-pisces.final.vcf.gz' -o -name '*-basespace-pisces.final.vcf' \\) | head -1)
    [ -n "\$VCF" ] || { echo "ERROR: No *-basespace-pisces.final.vcf[.gz] found in ${sample_folder}" >&2; exit 1; }

    cnv_depth_qc.py "\$VCF" "${meta.sample_id}" \\
        --min-depth ${params.cnv_min_depth} \\
        --min-fraction ${params.cnv_min_depth_fraction} \\
        ${pass_only} > "${meta.sample_id}_cnv_depth_qc.txt"
    """
}
