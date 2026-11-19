include { CLINICAL_PATIENTS    } from '../../../modules/local/clinical_patients/main'
include { CLINICAL_SAMPLES     } from '../../../modules/local/clinical_samples/main'
include { WRITE_CASE_LISTS     } from '../../../modules/local/write_case_lists/main'
include { WRITE_META           } from '../../../modules/local/write_meta/main'
include { WRITE_GENE_PANELS    } from '../../../modules/local/write_gene_panels/main'
include { ANONYMIZE_CLINICAL as ANON_PATIENT } from '../../../modules/local/anonymize_clinical/main'
include { ANONYMIZE_CLINICAL as ANON_SAMPLE  } from '../../../modules/local/anonymize_clinical/main'

workflow STUDY_METADATA {

    take:
    ch_patient_file   // channel: patient file path
    ch_sample_file    // channel: sample file path
    ch_orig_linking   // value channel: original filtered linking file (real IDs)
    ch_output_linking // value channel: linking file for output IDs (anon or orig)
    ch_cnv_failed     // channel: list of sample IDs that failed CNV depth QC
    study_id          // val: study ID string

    main:
    CLINICAL_PATIENTS(ch_patient_file, ch_orig_linking)
    CLINICAL_SAMPLES(ch_sample_file, ch_orig_linking)

    if (params.anonymize?.toString()?.toBoolean()) {
        ANON_PATIENT(CLINICAL_PATIENTS.out, ch_orig_linking, ch_output_linking)
        ANON_SAMPLE(CLINICAL_SAMPLES.out, ch_orig_linking, ch_output_linking)
        ch_clinical_patient = ANON_PATIENT.out
        ch_clinical_sample  = ANON_SAMPLE.out
    } else {
        ch_clinical_patient = CLINICAL_PATIENTS.out
        ch_clinical_sample  = CLINICAL_SAMPLES.out
    }

    WRITE_CASE_LISTS(ch_output_linking, ch_cnv_failed, study_id)
    WRITE_META(study_id)

    // Gene panel matrix: every sample profiled with the panel, so frequencies use the whole cohort
    ch_meta_files  = WRITE_META.out
    ch_panel_data  = Channel.empty()
    ch_gene_panels = Channel.value([])
    if (params.gene_panel && params.gene_panel.toString() != 'null') {
        WRITE_GENE_PANELS(file(params.gene_panel, checkIfExists: true), ch_output_linking, ch_cnv_failed, study_id)
        ch_meta_files  = WRITE_META.out.mix(WRITE_GENE_PANELS.out.meta).collect()
        ch_panel_data  = WRITE_GENE_PANELS.out.matrix
        ch_gene_panels = WRITE_GENE_PANELS.out.panels
    }

    emit:
    clinical_patient = ch_clinical_patient
    clinical_sample  = ch_clinical_sample
    meta_files       = ch_meta_files
    case_lists       = WRITE_CASE_LISTS.out
    gene_panel_data  = ch_panel_data
    gene_panels      = ch_gene_panels
}
