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
    if [ -n "\$VCF" ]; then
        vcf_to_seg.py "\$VCF" "${meta.sample_id}"
    else
        printf 'ID\\tchrom\\tloc.start\\tloc.end\\tnum.mark\\tseg.mean\\n' > data_seg.txt
    fi
    mv data_seg.txt "${meta.sample_id}_seg.txt"
    """
}
