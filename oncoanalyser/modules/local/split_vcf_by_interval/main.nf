process SPLIT_VCF_BY_INTERVAL {
    tag "${meta.sample}.${interval_meta.id}"

    label "process_medium"

    // take apptainer container from nf-core https://github.com/nf-core/modules/blob/master/modules/nf-core/bcftools/view/main.nf
    container 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/47/474a5ea8dc03366b04df884d89aeacc4f8e6d1ad92266888e7a8e7958d07cde8/data'

    input:
    tuple val(meta), path(vcf), path(index), val(interval_meta), path(interval_bed)

    output:
    tuple val(meta), val(interval_meta.id), path("${meta.id}.${interval_meta.id}.vcf.gz"), path("${meta.id}.${interval_meta.id}.vcf.gz.tbi")

    script:
    """
    # mv needed to match with index file
    mv ${vcf} ${meta.sample}.vcf.gz

    bcftools view -R ${interval_bed} ${meta.sample}.vcf.gz -O z -o ${meta.sample}.${interval_meta.id}.vcf.gz
    bcftools index -t ${meta.sample}.${interval_meta.id}.vcf.gz
    """
}
