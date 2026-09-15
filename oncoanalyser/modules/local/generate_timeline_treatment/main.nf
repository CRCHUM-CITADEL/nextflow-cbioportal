process GENERATE_TIMELINE_TREATMENT {
    publishDir { "${params.outdir}/${meta.group}" }, mode: 'copy', enabled: false

    container params.container_r

    tag { meta.group }

    input:
        tuple val(meta),
              path(sample_registrations),
              path(treatments),
              path(systemic_therapies),
              path(genomic_subjects),
              path(treatment_intent_map)
        path(clinical_common_r)

    output:
        tuple val(meta.group), path("data_timeline_treatment.txt"), emit: ch_timeline_part

    script:
    def sample_reg_arg       = sample_registrations ? "--sample_registrations ${sample_registrations}" : ""
    def treatments_arg       = treatments            ? "--treatments ${treatments}"                     : ""
    def systemic_arg         = systemic_therapies    ? "--systemic_therapies ${systemic_therapies}"     : ""
    def genomic_subjects_arg = genomic_subjects      ? "--genomic_subjects ${genomic_subjects}"         : ""
    """
    gen_timeline_treatment.R \
        ${sample_reg_arg} \
        ${treatments_arg} \
        ${systemic_arg} \
        ${genomic_subjects_arg} \
        --treatment_intent_map ${treatment_intent_map}
    """

    stub:
    """
    touch data_timeline_treatment.txt
    """
}
