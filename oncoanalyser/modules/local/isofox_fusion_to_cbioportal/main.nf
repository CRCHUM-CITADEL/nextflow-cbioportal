process ISOFOX_FUSION_TO_CBIOPORTAL {
    tag "$meta.sample"
    label 'process_single'

    container params.container_r

    input:
        tuple val(meta), path(fusion_csv)

    output:
        tuple val(meta), path("${meta.sample}.isofox_fusion.data_sv.txt"), emit: sv

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/gen_isofox_fusion_to_cbioportal.R \\
        --input  ${fusion_csv} \\
        --sample ${meta.sample} \\
        --output ${meta.sample}.isofox_fusion.data_sv.txt
    """

    stub:
    """
    touch ${meta.sample}.isofox_fusion.data_sv.txt
    """
}
