process GENERATE_CASE_LIST {
    publishDir "${params.outdir}/case_lists/", mode: 'copy'

    input:
    val label
    val(list_of_samples)

    output:
    path "cases_${label}.txt"

    script:
    """
    echo -e "cancer_study_identifier: add_text" > cases_${label}.txt
    echo -e "stable_id: add_text_${label}" >> cases_${label}.txt
    echo -e "case_list_name: add_text" >> cases_${label}.txt
    echo -e "case_list_description: ADD TEXT" >> cases_${label}.txt
    echo -e "case_list_ids: ${list_of_samples}" >> cases_${label}.txt
    """
}
