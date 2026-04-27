process SIGS_TO_CBIOPORTAL {
    tag "$meta.sample"
    label 'process_single'

    container params.container_r

    publishDir "${params.outdir}/${meta.group}/${meta.subject}", mode: 'copy'

    input:
        tuple val(meta), path(sig_alloc)

    output:
        tuple val(meta), path("${meta.sample}.data_sigs.txt"), emit: sigs

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/gen_sigs_to_cbioportal.R \\
        --input  ${sig_alloc} \\
        --sample ${meta.sample} \\
        --output ${meta.sample}.data_sigs.txt
    """

    stub:
    """
    touch ${meta.sample}.data_sigs.txt
    """
}
