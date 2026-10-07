process PACKAGE_CBIOPORTAL {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    val study_id
    path data_files
    path meta_files
    path case_lists

    output:
    path "${study_id}.tar.gz"

    script:
    """
    mkdir -p "${study_id}"
    cp ${data_files} "${study_id}/"
    cp ${meta_files} "${study_id}/"
    cp -rL case_lists "${study_id}/"   # -L: the staged case_lists is a symlink; copy its contents
    tar -czf "${study_id}.tar.gz" "${study_id}"
    """

    stub:
    """
    touch "${study_id}.tar.gz"
    """
}
