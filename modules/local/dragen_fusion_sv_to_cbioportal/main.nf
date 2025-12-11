process DRAGEN_FUSION_SV_TO_CBIOPORTAL {
    publishDir "${params.outdir}/${meta.group}/${meta.subject}", mode: 'copy'

    tag { meta.subject }

    container params.container_r

    input:
        tuple val(meta), path(dragen_fusion)

    output:
        tuple val(meta), path("${meta.subject}.data_sv.txt")

    script:
    """
    gen_format_dragen_fusion.R \
        -i $dragen_fusion \
        -o ${meta.subject}.data_sv.txt \
        -s ${meta.subject}
    """
}
