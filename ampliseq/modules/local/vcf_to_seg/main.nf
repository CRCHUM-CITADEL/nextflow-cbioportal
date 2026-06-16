process VCF_TO_SEG {
    tag "${meta.sample_id}"
    label 'python'
    publishDir "${params.outdir}/samples/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(sample_folder)

    output:
    path("${meta.sample_id}_seg.txt")

    script:
    """
    VCF=\$(find -L "${sample_folder}" -maxdepth 1 -name '*-basespace-cnv.final.vcf' | head -1)
    [ -n "\$VCF" ] || { echo "ERROR: No *-basespace-cnv.final.vcf found in ${sample_folder}" >&2; exit 1; }

    vcf_to_seg.py "\$VCF" "${meta.subject_id}"
    mv data_seg.txt "${meta.sample_id}_seg.txt"
    """
}
