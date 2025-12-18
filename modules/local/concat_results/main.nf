process CONCAT_RESULTS {
    input:
    tuple val(meta), path(vcfs), path(indexes)
    
    output:
    tuple val(meta), path("${meta.id}.merged.vcf.gz"), path("${meta.id}.merged.vcf.gz.tbi")
    
    script:
    """
    bcftools concat ${vcfs} -O z -o ${meta.id}.merged.vcf.gz
    bcftools index -t ${meta.id}.merged.vcf.gz
    """
}