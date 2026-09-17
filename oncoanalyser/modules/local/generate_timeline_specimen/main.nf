process GENERATE_TIMELINE_SPECIMEN {
    publishDir { "${params.outdir}/${meta.group}" }, mode: 'copy', enabled: false

    container params.container_r

    tag { meta.group }

    input:
        tuple val(meta),
              path(sample_registrations),
              path(specimens),
              path(follow_ups),
              path(genomic_subjects),
              path(primary_site_map),
              path(specimen_tissue_source_map)
        path(clinical_common_r)

    output:
        tuple val(meta.group), path("data_timeline_specimen.txt"), emit: ch_timeline_part, optional: true

    script:
    def sample_reg_arg       = sample_registrations ? "--sample_registrations ${sample_registrations}" : ""
    def specimens_arg        = specimens             ? "--specimens ${specimens}"                       : ""
    def follow_ups_arg       = follow_ups            ? "--follow_ups ${follow_ups}"                     : ""
    def genomic_subjects_arg = genomic_subjects      ? "--genomic_subjects ${genomic_subjects}"         : ""
    """
    gen_timeline_specimen.R \
        ${sample_reg_arg} \
        ${specimens_arg} \
        ${follow_ups_arg} \
        ${genomic_subjects_arg} \
        --primary_site_map ${primary_site_map} \
        --specimen_tissue_source_map ${specimen_tissue_source_map}
    """

    stub:
    """
    touch data_timeline_specimen.txt
    """
}
