include { CLINICAL_PATIENTS } from '../../../modules/local/clinical_patients/main'
include { CLINICAL_SAMPLES  } from '../../../modules/local/clinical_samples/main'
include { WRITE_CASE_LISTS  } from '../../../modules/local/write_case_lists/main'
include { WRITE_META        } from '../../../modules/local/write_meta/main'

workflow STUDY_METADATA {

    take:
    ch_patient_file   // channel: patient file path
    ch_sample_file    // channel: sample file path
    ch_linking        // value channel: filtered linking file (samplesheet samples only)
    study_id          // val: study ID string

    main:
    CLINICAL_PATIENTS(ch_patient_file, ch_linking)
    CLINICAL_SAMPLES(ch_sample_file, ch_linking)
    WRITE_CASE_LISTS(ch_linking, study_id)
    WRITE_META(study_id)

    emit:
    clinical_patient = CLINICAL_PATIENTS.out
    clinical_sample  = CLINICAL_SAMPLES.out
    meta_files       = WRITE_META.out
    case_lists       = WRITE_CASE_LISTS.out
}
