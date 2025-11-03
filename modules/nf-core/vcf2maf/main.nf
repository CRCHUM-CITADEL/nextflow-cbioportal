// modified to allow vcf.gz and vcf (10/10/2025)
process VCF2MAF {
    tag "$meta.id"
    label 'process_medium_memory'
    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    conda "${moduleDir}/environment.yml"
    // added local container
    container "${ workflow.containerEngine == 'apptainer' && !task.ext.singularity_pull_docker_container ?
        (params.container_vcf2maf ?: 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/7c/7cbf9421f0bee23a93a35c5d0c7166ac1e89a40008d8e474cecfddb93226bf65/data') :
        'community.wave.seqera.io/library/ensembl-vep_vcf2maf:2d40b60b4834af73' }"

    input:
        tuple val(meta), path(vcf) // Now accepts both compressed (.vcf.gz) and uncompressed (.vcf) files
        path fasta                 // Required
        path vep_cache             // Required for VEP running. A default of /.vep is supplied.

    output:
        tuple val(meta), path("*.maf"), emit: maf
        path "versions.yml"           , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args   ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    def vep_cache_cmd = vep_cache       ? "--vep-data $vep_cache ${params.vep_params}" : ""     // If VEP is present, it will find it and add it to commands otherwise blank
    def VERSION       = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ "$vep_cache" ]; then
        VEP_CMD="--vep-path \$(dirname \$(type -p vep))"
        VEP_VERSION=\$(echo -e "\\n    ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')")
    else
        VEP_CMD=""
        VEP_VERSION=""
    fi

    echo -e "\$VEP_VERSION"

    # Handle compressed VCF files
    if [[ $vcf == *.gz ]]; then
        tmp=\$(mktemp --suffix=.vcf)
        rm -f "\$tmp" 
        gunzip -c "$vcf" > "\$tmp"
        INPUT_VCF="\$tmp"
    else
        INPUT_VCF="$vcf"
    fi

    originalTumorName=`cat \$INPUT_VCF| grep CHROM | awk -F"\\t" '{print \$(NF)}'`
    originalNormalName=`cat \$INPUT_VCF | grep CHROM | awk -F"\\t" '{print \$(NF-1)}'`

    vcf2maf.pl \\
        $args \\
        --tumor-id \$originalTumorName \\
        --normal-id \$originalNormalName \\
        \$VEP_CMD \\
        $vep_cache_cmd \\
        --ref-fasta $fasta \\
        --input-vcf \$INPUT_VCF \\
        --output-maf ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION\$VEP_VERSION
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def VERSION = '1.6.22' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    if [ "$vep_cache" ]; then
        VEP_VERSION=\$(echo -e "\\n    ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')")
    else
        VEP_VERSION=""
    fi

    touch ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: $VERSION\$VEP_VERSION
    END_VERSIONS
    """
}
