process GENERATE_TIMELINE_STATUS {
    publishDir { "${params.outdir}/${meta.group}" }, mode: 'copy', enabled: false

    container params.container_r

    tag { meta.group }

    input:
        tuple val(meta),
              path(sample_registrations),
              path(follow_ups),
              path(genomic_subjects),
              path(primary_site_map)
        path(clinical_common_r)

    output:
        tuple val(meta.group), path("data_timeline_status.txt"), emit: ch_timeline_part

    script:
    def sample_reg_arg       = sample_registrations ? "--sample_registrations ${sample_registrations}" : ""
    def follow_ups_arg       = follow_ups            ? "--follow_ups ${follow_ups}"                     : ""
    def genomic_subjects_arg = genomic_subjects      ? "--genomic_subjects ${genomic_subjects}"         : ""
    """
    gen_timeline_status.R \
        ${sample_reg_arg} \
        ${follow_ups_arg} \
        ${genomic_subjects_arg} \
        --primary_site_map ${primary_site_map}
    """

    stub:
    """
    touch data_timeline_status.txt
    """
}
