process FORMAT_SV {
    tag "${meta.sample_id}"
    label 'python'
    publishDir { "${params.outdir}/samples/${meta.sample_id}" }, mode: 'copy'

    input:
    tuple val(meta), path(tsv), path(sample_folder)
    path(gene_loci) // GRCh37 gene spans, places the TSV fusion partner's chromosome

    output:
    path("${meta.sample_id}_sv.txt")

    script:
    // Filter mode: always the TSV (its Supporting Reads), even with a STAR-Fusion VCF
    if (params.filter_tsv_variants.toString().toBoolean()) {
        """
        format_tsv.py "${tsv}" "${meta.sample_id}" --gene-loci "${gene_loci}" --min-supporting-reads ${params.sv_min_supporting_reads}
        mv data_sv.txt "${meta.sample_id}_sv.txt"
        """
    } else {
        """
        VCF=\$(find -L "${sample_folder}" -maxdepth 1 -name '*-star-fusion.final.vcf' | head -1)
        if [ -n "\$VCF" ]; then
            fusion_vcf_to_sv.py "\$VCF" "${meta.sample_id}"
        else
            format_tsv.py "${tsv}" "${meta.sample_id}" --gene-loci "${gene_loci}"
        fi
        mv data_sv.txt "${meta.sample_id}_sv.txt"
        """
    }
}
