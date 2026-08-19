process FORMAT_CLINICAL {
    publishDir { "${params.outdir}/${meta.group}" }, mode: 'copy'

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
              path(genomic_subjects),
              path(primary_site_map),
              path(specimen_tissue_source_map),
              path(treatment_intent_map)

    output:
        tuple val(meta), path("data_clinical_${meta.mode}.txt")

    script:
    def treatments_arg        = treatments         ? "--treatments ${treatments}"                      : ""
    def surgeries_arg         = surgeries          ? "--surgeries ${surgeries}"                        : ""
    def systemic_arg          = systemic_therapies ? "--systemic_therapies ${systemic_therapies}"      : ""
    def radiations_arg        = radiations         ? "--radiations ${radiations}"                      : ""
    def follow_ups_arg        = follow_ups         ? "--follow_ups ${follow_ups}"                      : ""
    def biomarkers_arg        = biomarkers         ? "--biomarkers ${biomarkers}"                      : ""
    def genomic_subjects_arg      = genomic_subjects         ? "--genomic_subjects ${genomic_subjects}"                     : ""
    def primary_site_map_arg      = "--primary_site_map ${primary_site_map}"
    def specimen_source_map_arg   = "--specimen_tissue_source_map ${specimen_tissue_source_map}"
    def treatment_intent_map_arg  = "--treatment_intent_map ${treatment_intent_map}"
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
        ${primary_site_map_arg} \
        ${specimen_source_map_arg} \
        ${treatment_intent_map_arg} \
        --output data_clinical_${meta.mode}.txt
    """

    stub:
    """
    touch data_clinical_${meta.mode}.txt
    """
}
