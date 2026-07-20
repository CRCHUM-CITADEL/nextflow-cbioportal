process FORMAT_CLINICAL {
    publishDir "${params.outdir}/${meta.group}", mode: 'copy'

    container params.container_r

    tag { meta.mode + meta.group }

    input:
        tuple val(meta), val(sample_list)

    output:
        path "data_clinical_${meta.mode}.txt"

    script:
    def cancer_type           = params.oncotree_code ? params.oncotree_code.toLowerCase() : ""
    def icd_arg               = params.icd_code ? "--icd_code ${params.icd_code}" : ""
    def treatments_arg        = sample_list.treatments        ? "--treatments ${sample_list.treatments}"               : ""
    def surgeries_arg         = sample_list.surgeries         ? "--surgeries ${sample_list.surgeries}"                 : ""
    def systemic_arg          = sample_list.systemic_therapies ? "--systemic_therapies ${sample_list.systemic_therapies}" : ""
    def radiations_arg        = sample_list.radiations        ? "--radiations ${sample_list.radiations}"               : ""
    def follow_ups_arg        = sample_list.follow_ups        ? "--follow_ups ${sample_list.follow_ups}"               : ""
    def biomarkers_arg        = sample_list.biomarkers        ? "--biomarkers ${sample_list.biomarkers}"               : ""
    """
    clin_format.R \
        --mode ${meta.mode} \
        --donors ${sample_list.donors} \
        --primary_diagnoses ${sample_list.primary_diagnoses} \
        --specimens ${sample_list.specimens} \
        --sample_registrations ${sample_list.sample_registrations} \
        ${treatments_arg} \
        ${surgeries_arg} \
        ${systemic_arg} \
        ${radiations_arg} \
        ${follow_ups_arg} \
        ${biomarkers_arg} \
        --cancer_type ${cancer_type} \
        ${icd_arg} \
        --output data_clinical_${meta.mode}.txt
    """
}
