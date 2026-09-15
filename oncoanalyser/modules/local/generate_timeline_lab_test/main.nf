process GENERATE_TIMELINE_LAB_TEST {
    publishDir { "${params.outdir}/${meta.group}" }, mode: 'copy', enabled: false

    container params.container_r

    tag { meta.group }

    input:
        tuple val(meta),
              path(sample_registrations),
              path(biomarkers),
              path(genomic_subjects)
        path(clinical_common_r)

    output:
        tuple val(meta.group), path("data_timeline_lab_test.txt"), emit: ch_timeline_part

    script:
    def sample_reg_arg       = sample_registrations ? "--sample_registrations ${sample_registrations}" : ""
    def biomarkers_arg       = biomarkers            ? "--biomarkers ${biomarkers}"                     : ""
    def genomic_subjects_arg = genomic_subjects      ? "--genomic_subjects ${genomic_subjects}"         : ""
    """
    gen_timeline_lab_test.R \
        ${sample_reg_arg} \
        ${biomarkers_arg} \
        ${genomic_subjects_arg}
    """

    stub:
    """
    touch data_timeline_lab_test.txt
    """
}
