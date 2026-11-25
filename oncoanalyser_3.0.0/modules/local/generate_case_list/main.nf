process GENERATE_CASE_LIST {
    publishDir { "${params.outdir}/${group}/case_lists" }, mode: 'copy'

    input:
    tuple val(group), val(suffix), val(name), val(description), val(category), val(samples)

    output:
    tuple val(group), path("cases_${suffix}.txt")

    script:
    def group_lower = group.toLowerCase().replace('-', '_')
    // quoted heredoc: the shell expands nothing, so the tabs between sample IDs survive
    """
    cat <<'EOF' > cases_${suffix}.txt
cancer_study_identifier: ${group_lower}
stable_id: ${group_lower}_${suffix}
case_list_name: ${name}
case_list_description: ${description}
case_list_category: ${category}
case_list_ids: ${samples}
EOF
    """

    stub:
    """
    touch "cases_${suffix}.txt"
    """
}
