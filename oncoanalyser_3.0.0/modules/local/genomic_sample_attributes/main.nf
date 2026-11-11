process GENOMIC_SAMPLE_ATTRIBUTES {
    tag "${meta.sample}"
    label 'process_single'

    container "${params.container_python}"

    input:
        tuple val(meta), path(purity), path(chord) // [] when absent

    output:
        tuple val(meta), path("${meta.sample}.genomic_sample_attributes.tsv"), emit: attributes

    when:
        task.ext.when == null || task.ext.when

    script:
    def purity_arg = purity ? "--purity ${purity}" : ''
    def chord_arg  = chord  ? "--chord ${chord}"   : ''
    """
    gen_genomic_sample_attributes.py \\
        --patient ${meta.subject} \\
        --sample ${meta.sample} \\
        ${purity_arg} \\
        ${chord_arg} \\
        --output ${meta.sample}.genomic_sample_attributes.tsv
    """

    stub:
    """
    {
        printf '#Patient Identifier\\tSample Identifier\\tTMB per Mb (PURPLE)\\tTMB Status (PURPLE)\\tTML (PURPLE)\\tTML Status (PURPLE)\\tHRD Score (CHORD)\\tHR Status (CHORD)\\tHRD Type (CHORD)\\n'
        printf '#stub\\tstub\\tstub\\tstub\\tstub\\tstub\\tstub\\tstub\\tstub\\n'
        printf '#STRING\\tSTRING\\tNUMBER\\tSTRING\\tNUMBER\\tSTRING\\tNUMBER\\tSTRING\\tSTRING\\n'
        printf '#1\\t1\\t1\\t1\\t1\\t1\\t1\\t1\\t1\\n'
        printf 'PATIENT_ID\\tSAMPLE_ID\\tTMB_PER_MB\\tTMB_STATUS\\tTML\\tTML_STATUS\\tHRD_SCORE\\tHR_STATUS\\tHRD_TYPE\\n'
        printf '%s\\t%s\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA\\n' "${meta.subject}" "${meta.sample}"
    } > ${meta.sample}.genomic_sample_attributes.tsv
    """
}
