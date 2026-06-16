include { CLINICAL_PATIENTS    } from '../../../modules/local/clinical_patients/main'
include { CLINICAL_SAMPLES     } from '../../../modules/local/clinical_samples/main'
include { WRITE_CASE_LISTS     } from '../../../modules/local/write_case_lists/main'
include { WRITE_META           } from '../../../modules/local/write_meta/main'
include { ANONYMIZE_CLINICAL as ANON_PATIENT } from '../../../modules/local/anonymize_clinical/main'
include { ANONYMIZE_CLINICAL as ANON_SAMPLE  } from '../../../modules/local/anonymize_clinical/main'

workflow STUDY_METADATA {

    take:
    ch_patient_file   // channel: patient file path
    ch_sample_file    // channel: sample file path
    ch_linking        // value channel: filtered linking file (samplesheet samples only)
    ch_patient_map    // value channel: subject_id<TAB>sample_id mapping file
    study_id          // val: study ID string

    main:
    CLINICAL_PATIENTS(ch_patient_file, ch_linking)
    CLINICAL_SAMPLES(ch_sample_file, ch_linking)

    if (params.anonymize) {
        ANON_PATIENT(CLINICAL_PATIENTS.out, ch_linking, ch_patient_map)
        ANON_SAMPLE(CLINICAL_SAMPLES.out, ch_linking, ch_patient_map)
        ch_clinical_patient = ANON_PATIENT.out
        ch_clinical_sample  = ANON_SAMPLE.out
    } else {
        ch_clinical_patient = CLINICAL_PATIENTS.out
        ch_clinical_sample  = CLINICAL_SAMPLES.out
    }

    WRITE_CASE_LISTS(ch_linking, study_id)
    WRITE_META(study_id)

    emit:
    clinical_patient = ch_clinical_patient
    clinical_sample  = ch_clinical_sample
    meta_files       = WRITE_META.out
    case_lists       = WRITE_CASE_LISTS.out
}
