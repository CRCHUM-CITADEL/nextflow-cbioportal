process ESVEE_SV_TO_CBIOPORTAL {
    tag "$meta.sample"
    label 'process_low'

    container params.container_r

    publishDir "${params.outdir}/${meta.group}/${meta.subject}", mode: 'copy'

    input:
        tuple val(meta), path(esvee_vcf_tumor)
        path ensembl_annotations

    output:
        tuple val(meta), path("${meta.sample}.data_sv.txt"), emit: sv

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/gen_esvee_sv_to_cbioportal.R \\
        --input               ${esvee_vcf_tumor} \\
        --sample              ${meta.sample} \\
        --ensembl_annotations ${ensembl_annotations} \\
        --output              ${meta.sample}.data_sv.txt
    """

    stub:
    """
    touch ${meta.sample}.data_sv.txt
    """
}
