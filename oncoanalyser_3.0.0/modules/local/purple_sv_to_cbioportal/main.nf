process PURPLE_SV_TO_CBIOPORTAL {
    tag "$meta.sample"
    label 'process_low'

    container params.container_r

    input:
        tuple val(meta), path(purple_sv_vcf)
        path ensembl_annotations

    output:
        tuple val(meta), path("${meta.sample}.purple_sv.data_sv.txt"), emit: sv

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    gen_purple_sv_to_cbioportal.R \\
        --input               ${purple_sv_vcf} \\
        --sample              ${meta.sample} \\
        --tumor_id            ${meta.subject}-T \\
        --ensembl_annotations ${ensembl_annotations} \\
        --output              ${meta.sample}.purple_sv.data_sv.txt
    """

    stub:
    """
    touch ${meta.sample}.purple_sv.data_sv.txt
    """
}
