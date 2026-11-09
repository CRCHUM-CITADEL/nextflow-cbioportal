process PACKAGE_CBIOPORTAL {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val study_id
    path data_files
    path meta_files
    path case_lists
    path gene_panels   // panel definitions for importGenePanel.pl, or [] without gene_panel

    output:
    path "${study_id}.tar.gz"

    script:
    """
    mkdir -p "${study_id}"
    cp ${data_files} "${study_id}/"
    cp ${meta_files} "${study_id}/"
    cp -rL case_lists "${study_id}/"   # -L: the staged case_lists is a symlink; copy its contents
    if [ -n "${gene_panels}" ]; then cp -rL ${gene_panels} "${study_id}/"; fi
    tar -czf "${study_id}.tar.gz" "${study_id}"
    """

    stub:
    """
    touch "${study_id}.tar.gz"
    """
}
