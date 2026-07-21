process FORMAT_CLINICAL {
    publishDir "${params.outdir}/${meta.group}", mode: 'copy'

    container params.container_r

    tag { meta.mode + meta.group }

    input:
        tuple val(meta),
              path(donors),
              path(primary_diagnoses),
              path(specimens),
              path(sample_registrations),
              path(treatments),
              path(surgeries),
              path(systemic_therapies),
              path(radiations),
              path(follow_ups),
              path(biomarkers),
              path(genomic_subjects)

    output:
        path "data_clinical_${meta.mode}.txt"

    script:
    def treatments_arg        = treatments         ? "--treatments ${treatments}"                      : ""
    def surgeries_arg         = surgeries          ? "--surgeries ${surgeries}"                        : ""
    def systemic_arg          = systemic_therapies ? "--systemic_therapies ${systemic_therapies}"      : ""
    def radiations_arg        = radiations         ? "--radiations ${radiations}"                      : ""
    def follow_ups_arg        = follow_ups         ? "--follow_ups ${follow_ups}"                      : ""
    def biomarkers_arg        = biomarkers         ? "--biomarkers ${biomarkers}"                      : ""
    def genomic_subjects_arg  = genomic_subjects   ? "--genomic_subjects ${genomic_subjects}"          : ""
    """
    clin_format.R \
        --mode ${meta.mode} \
        --donors ${donors} \
        --primary_diagnoses ${primary_diagnoses} \
        --specimens ${specimens} \
        --sample_registrations ${sample_registrations} \
        ${treatments_arg} \
        ${surgeries_arg} \
        ${systemic_arg} \
        ${radiations_arg} \
        ${follow_ups_arg} \
        ${biomarkers_arg} \
        ${genomic_subjects_arg} \
        --output data_clinical_${meta.mode}.txt
    """

    stub:
    """
    touch data_clinical_${meta.mode}.txt
    """
}
