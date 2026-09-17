process PACKAGE_CBIOPORTAL {
    container params.container_r
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(group), path(all_files)

    output:
    path "${group}.tar.gz"

    script:
    """
    mkdir -p "${group}/case_lists"
    for f in ${all_files}; do
        if [[ \$(basename "\$f") == cases_* ]]; then
            cp "\$f" "${group}/case_lists/"
        else
            cp "\$f" "${group}/"
        fi
    done
    tar -czf "${group}.tar.gz" "${group}"
    """

    stub:
    """
    touch "${group}.tar.gz"
    """
}
