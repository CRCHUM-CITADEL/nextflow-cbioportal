process WRITE_CASE_LISTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(linking_file)
    val(study_id)

    output:
    path("case_lists")

    script:
    """
    mkdir -p case_lists

    IDS_TAB=\$(awk 'NR>1 {printf "%s\\t", \$2}' ${linking_file} | sed 's/\\t\$//')

    {
        echo "cancer_study_identifier: ${study_id}"
        echo "stable_id: ${study_id}_sequenced"
        echo "case_list_name: all mutations"
        echo "case_list_description: all mutations of ${study_id}"
        echo "case_list_ids: \${IDS_TAB}"
    } > case_lists/cases_sequenced.txt

    {
        echo "cancer_study_identifier: ${study_id}"
        echo "stable_id: ${study_id}_sv"
        echo "case_list_name: all sv"
        echo "case_list_description: all sv of ${study_id}"
        echo "case_list_ids: \${IDS_TAB}"
    } > case_lists/cases_sv.txt

    {
        echo "cancer_study_identifier: ${study_id}"
        echo "stable_id: ${study_id}_cna"
        echo "case_list_name: all cna"
        echo "case_list_description: all cna of ${study_id}"
        echo "case_list_ids: \${IDS_TAB}"
    } > case_lists/cases_cna.txt
    """
}
