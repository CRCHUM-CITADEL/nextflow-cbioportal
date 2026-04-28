process MERGE_SIGS_TO_CBIOPORTAL {
    publishDir "${params.outdir}/${meta.group}", mode: 'copy'

    container params.container_r

    input:
        tuple val(meta), path(sigs_file_list)

    output:
        tuple val(meta.group), path("data_mutational_signatures_contribution_SBS.txt")

    when:
        task.ext.when == null || task.ext.when

    script:
    """
    Rscript ${projectDir}/bin/gen_merge_sigs_to_cbioportal.R \\
        --input_files ${sigs_file_list.join(',')} \\
        --output_file data_mutational_signatures_contribution_SBS.txt
    """

    stub:
    """
    touch data_mutational_signatures_contribution_SBS.txt
    """
}
