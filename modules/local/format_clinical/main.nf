process FORMAT_CLINICAL {
    publishDir "${params.outdir}/", mode: 'copy'

    container params.container_r

    tag { meta.mode + meta.group }

    input:
        tuple val(meta), val(sample_list)

    output:
        path "data_clinical_${meta.mode}.txt"

    script:
    """
    clin_format.R \
        --mode ${meta.mode} \
        --patient ${sample_list.patient} \
        --diagnosis ${sample_list.diagnosis} \
        --treatment ${sample_list.treatment} \
        --surgeries ${sample_list.surgeries} \
        --systreat ${sample_list.systemic_treatment} \
        --specimen ${sample_list.specimen} \
        --radiotherapy ${sample_list.radio_therapy} \
        --output data_clinical_${meta.mode}.txt
    """
}
