process FORMAT_CLINICAL {
    publishDir "${params.outdir}/${meta.group}", mode: 'copy'

    container params.container_r

    tag { meta.mode + meta.group }

    input:
        tuple val(meta), val(sample_list)

    output:
        path "data_clinical_${meta.mode}.txt"

    script:
    def treatments_arg        = sample_list.treatments        ? "--treatments ${sample_list.treatments}"               : ""
    def surgeries_arg         = sample_list.surgeries         ? "--surgeries ${sample_list.surgeries}"                 : ""
    def systemic_arg          = sample_list.systemic_therapies ? "--systemic_therapies ${sample_list.systemic_therapies}" : ""
    def radiations_arg        = sample_list.radiations        ? "--radiations ${sample_list.radiations}"               : ""
    def follow_ups_arg        = sample_list.follow_ups        ? "--follow_ups ${sample_list.follow_ups}"               : ""
    def biomarkers_arg        = sample_list.biomarkers        ? "--biomarkers ${sample_list.biomarkers}"               : ""
    def genomic_subjects_arg  = sample_list.genomic_subjects  ? "--genomic_subjects ${sample_list.genomic_subjects}"  : ""
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
        ${genomic_subjects_arg} \
        --output data_clinical_${meta.mode}.txt
    """
}
