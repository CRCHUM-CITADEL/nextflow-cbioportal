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
    """
    zcat $ger_dna_vcf | grep "#" > tmp.${meta.sample}.vcf
    zcat $ger_dna_vcf | grep PASS | grep -v "FILTER=<ID=low_depth" >> tmp.${meta.sample}.vcf

    bcftools view \\
        -s ^${meta.sample} \\
        -Oz \\
        -o ${meta.sample}.vcf.gz \\
        tmp.${meta.sample}.vcf
    """
}
