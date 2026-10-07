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

    # case_list_category marks these as cBioPortal's standard lists (OncoPrint and the query
    # page pick them by category); every samplesheet sample is profiled for all three types.
    write_list() {  # <file> <stable_id suffix> <category> <name> <description>
        {
            echo "cancer_study_identifier: ${study_id}"
            echo "stable_id: ${study_id}_\$2"
            echo "case_list_category: \$3"
            echo "case_list_name: \$4"
            echo "case_list_description: \$5"
            echo "case_list_ids: \${IDS_TAB}"
        } > "case_lists/\$1"
    }
    write_list cases_sequenced.txt sequenced all_cases_with_mutation_data "Samples with mutation data" "All samples with mutation data in ${study_id}"
    write_list cases_cna.txt cna all_cases_with_cna_data "Samples with CNA data" "All samples with copy number alteration data in ${study_id}"
    write_list cases_sv.txt sv all_cases_with_sv_data "Samples with SV data" "All samples with structural variant data in ${study_id}"
    write_list cases_cnaseq.txt cnaseq all_cases_with_mutation_and_cna_data "Samples with mutation and CNA data" "All samples with both mutation and copy number alteration data in ${study_id}"
    """
}
