process FORMAT_SV {
    tag "${meta.sample_id}"
    label 'python'
    publishDir "${params.outdir}/samples/${meta.sample_id}", mode: 'copy'

    input:
    tuple val(meta), path(tsv), path(sample_folder)

    output:
    path("${meta.sample_id}_sv.txt")

    script:
    """
    VCF=\$(find -L "${sample_folder}" -maxdepth 1 -name '*-star-fusion.final.vcf' | head -1)
    if [ -n "\$VCF" ]; then
        fusion_vcf_to_sv.py "\$VCF" "${meta.sample_id}"
    else
        format_tsv.py "${tsv}" "${meta.sample_id}"
    fi
    mv data_sv.txt "${meta.sample_id}_sv.txt"
    """
}
