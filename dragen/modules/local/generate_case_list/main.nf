process GENERATE_CASE_LIST {
    publishDir "${params.outdir}/${group}/case_lists", mode: 'copy'

    input:
    val group
    val label
    val list_of_samples

    output:
    path "cases_${label}.txt"

    script:
    def group_lower = group.toLowerCase()
    """
    echo -e "cancer_study_identifier: ${group_lower}" > cases_${label}.txt
    echo -e "stable_id: ${group_lower}_${label}" >> cases_${label}.txt
    echo -e "case_list_name: add_text" >> cases_${label}.txt
    echo -e "case_list_description: ADD TEXT" >> cases_${label}.txt
    echo -e "case_list_ids: ${list_of_samples}" >> cases_${label}.txt
    """
}
