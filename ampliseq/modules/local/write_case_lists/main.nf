process WRITE_CASE_LISTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(linking_file)
    val(cnv_failed_ids) // sample IDs that failed CNV depth QC: left out of the CNA lists (may be empty)
    val(study_id)

    output:
    path("case_lists")

    script:
    """
    mkdir -p case_lists

    IDS_TAB=\$(awk 'NR>1 {printf "%s\\t", \$2}' ${linking_file} | sed 's/\\t\$//')
    CNA_IDS_TAB=\$(awk -v failed="${cnv_failed_ids.join(' ')}" '
        BEGIN { n = split(failed, ids, " "); for (i = 1; i <= n; i++) skip[ids[i]] = 1 }
        NR>1 && !(\$1 in skip) {printf "%s\\t", \$2}' ${linking_file} | sed 's/\\t\$//')

    # case_list_category marks these as cBioPortal's standard lists (OncoPrint and the query
    # page pick them by category). Every samplesheet sample is profiled for mutations and SVs;
    # samples that failed CNV depth QC are left out of the CNA lists (shown as not profiled).
    write_list() {  # <file> <stable_id suffix> <category> <name> <description> <ids>
        {
            echo "cancer_study_identifier: ${study_id}"
            echo "stable_id: ${study_id}_\$2"
            echo "case_list_category: \$3"
            echo "case_list_name: \$4"
            echo "case_list_description: \$5"
            echo "case_list_ids: \$6"
        } > "case_lists/\$1"
    }
    write_list cases_sequenced.txt sequenced all_cases_with_mutation_data "Samples with mutation data" "All samples with mutation data in ${study_id}" "\${IDS_TAB}"
    write_list cases_cna.txt cna all_cases_with_cna_data "Samples with CNA data" "All samples with copy number alteration data in ${study_id}" "\${CNA_IDS_TAB}"
    write_list cases_sv.txt sv all_cases_with_sv_data "Samples with SV data" "All samples with structural variant data in ${study_id}" "\${IDS_TAB}"
    write_list cases_cnaseq.txt cnaseq all_cases_with_mutation_and_cna_data "Samples with mutation and CNA data" "All samples with both mutation and copy number alteration data in ${study_id}" "\${CNA_IDS_TAB}"
    """
}
