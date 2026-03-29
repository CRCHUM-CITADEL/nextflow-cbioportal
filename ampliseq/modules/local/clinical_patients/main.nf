process CLINICAL_PATIENTS {
    label 'python'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(patient_file)

    output:
    path("data_clinical_patient.txt")

    script:
    """
    clinical_patients_format.py ${patient_file}
    """
}
