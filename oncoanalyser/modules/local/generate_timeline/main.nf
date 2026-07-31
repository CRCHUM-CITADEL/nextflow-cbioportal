process GENERATE_TIMELINE {
    publishDir { "${params.outdir}/${meta.group}" }, mode: 'copy'

    container params.container_r

    tag { meta.group }

    input:
        tuple val(meta),
              path(sample_registrations),
              path(treatments),
              path(surgeries),
              path(systemic_therapies),
              path(follow_ups),
              path(specimens),
              path(biomarkers),
              path(genomic_subjects),
              path(specimen_tissue_source_map),
              path(treatment_intent_map)

    output:
        tuple val(meta.group), path("data_timeline_*.txt", optional: true), emit: ch_timeline

    script:
    def sample_reg_arg            = sample_registrations    ? "--sample_registrations ${sample_registrations}"         : ""
    def treatments_arg            = treatments               ? "--treatments ${treatments}"                             : ""
    def surgeries_arg             = surgeries                ? "--surgeries ${surgeries}"                               : ""
    def systemic_arg              = systemic_therapies       ? "--systemic_therapies ${systemic_therapies}"             : ""
    def follow_ups_arg            = follow_ups               ? "--follow_ups ${follow_ups}"                             : ""
    def specimens_arg             = specimens                ? "--specimens ${specimens}"                                : ""
    def biomarkers_arg            = biomarkers               ? "--biomarkers ${biomarkers}"                             : ""
    def genomic_subjects_arg      = genomic_subjects         ? "--genomic_subjects ${genomic_subjects}"                 : ""
    def specimen_source_map_arg   = "--specimen_tissue_source_map ${specimen_tissue_source_map}"
    def treatment_intent_map_arg  = "--treatment_intent_map ${treatment_intent_map}"
    """
    gen_timeline.R \
        ${sample_reg_arg} \
        ${treatments_arg} \
        ${surgeries_arg} \
        ${systemic_arg} \
        ${follow_ups_arg} \
        ${specimens_arg} \
        ${biomarkers_arg} \
        ${genomic_subjects_arg} \
        ${specimen_source_map_arg} \
        ${treatment_intent_map_arg}
    """

    stub:
    """
    touch data_timeline_surgery.txt
    touch data_timeline_treatment.txt
    touch data_timeline_status.txt
    touch data_timeline_specimen.txt
    touch data_timeline_lab_test.txt
    """
}
