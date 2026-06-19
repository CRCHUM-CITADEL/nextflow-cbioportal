process VCF_TO_MAF {
    tag "${meta.sample_id}"
    label "vcf2maf"

    input:
    tuple val(meta), path(sample_folder)
    path(vep_data)
    path(ref_fasta)

    output:
    tuple val(meta), path("${meta.sample_id}.maf")

    script:
    """
    VCF_GZ=\$(find -L "${sample_folder}" -maxdepth 1 -name '*-basespace-pisces.final.vcf.gz' | head -1)
    [ -n "\$VCF_GZ" ] || { echo "ERROR: No VCF found in ${sample_folder}" >&2; exit 1; }

    gzip -dc "\$VCF_GZ" > "${meta.sample_id}.final.vcf"

    MAF_OUT="${meta.sample_id}-basespace-pisces.final.maf"

    vcf2maf.pl \\
        --tumor-id "${meta.sample_id}" --cache-version 113 \\
        --ncbi-build GRCh37 --vep-path /opt/conda/bin \\
        --vep-data ${vep_data} --ref-fasta ${ref_fasta} \\
        --input-vcf "${meta.sample_id}.final.vcf" \\
        --output-maf "\$MAF_OUT"

    # Strip vcf2maf version comment (first line); keep column header as line 1
    sed '1d' "\$MAF_OUT" > "${meta.sample_id}.maf"
    """
}
