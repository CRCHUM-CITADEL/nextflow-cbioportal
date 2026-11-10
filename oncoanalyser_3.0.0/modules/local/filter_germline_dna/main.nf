process FILTER_GERMLINE_DNA {
    tag "$meta.sample"
    label 'process_low'

    container "${ workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data':
        'community.wave.seqera.io/library/bcftools_htslib:0a3fa2654b52006f' }"

    input:
        tuple val(meta), path(ger_dna_vcf)

    output:
        tuple val(meta), path("*.vcf.gz")

    script:
    def min_dp  = params.germline_min_depth
    def het_vaf = params.germline_het_min_vaf
    def hom_vaf = params.germline_hom_min_vaf
    """
    # Keep FILTER exactly PASS or LOW_TUMOR_VCN (a tumour-side flag; the call in the normal is real), drop LOW_VAF and the rest.
    zcat $ger_dna_vcf | awk -F'\\t' '/^#/ || \$7 == "PASS" || \$7 == "LOW_TUMOR_VCN"' > tmp.${meta.sample}.vcf

    if ! grep -q '^##FORMAT=<ID=AF,' tmp.${meta.sample}.vcf || ! grep -q '^##FORMAT=<ID=DP,' tmp.${meta.sample}.vcf; then
        echo "ERROR: ${ger_dna_vcf} has no FORMAT/AF or FORMAT/DP, so germline genotypes cannot be inferred" >&2
        exit 1
    fi

    # PURPLE calls GT where it can and leaves ./. elsewhere; CPSR drops uncalled genotypes, so infer those from DP and AF.
    # Missing GTs only (PURPLE's calls are kept), in order: hom, het, then AF < het cutoff as reference (assumes no tumour-in-normal).
    # Low-depth sites with AF >= het cutoff stay ./. and CPSR drops them.
    # 0/0 is removed here: CPSR 2.2.5 keeps and classifies hom-ref records (only 2.3+ drops them).
    bcftools view -s ^${meta.sample} tmp.${meta.sample}.vcf \\
        | bcftools +setGT - -- -t q -n c:'1/1' -i 'GT="mis" & FMT/DP>=${min_dp} & FMT/AF>=${hom_vaf}' \\
        | bcftools +setGT - -- -t q -n c:'0/1' -i 'GT="mis" & FMT/DP>=${min_dp} & FMT/AF>=${het_vaf}' \\
        | bcftools +setGT - -- -t q -n c:'0/0' -i 'GT="mis" & FMT/AF<${het_vaf}' \\
        | bcftools view -e 'GT="ref"' -Oz -o ${meta.sample}.vcf.gz
    """

    stub:
    """
    echo | gzip > "${meta.sample}.vcf.gz"
    """
}
