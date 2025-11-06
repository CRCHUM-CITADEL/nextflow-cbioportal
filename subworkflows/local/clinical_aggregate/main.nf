include { FORMAT_CLINICAL } from '../../../modules/local/format_clinical'
include { ASSIGN_DATE } from '../../../modules/local/assign_date'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'

workflow CLINICAL_AGGREGATE {
    take:
        filelist
        id_linking_file

    main:

        csvs_with_date = ASSIGN_DATE(
                filelist
            ).map { meta, csv ->
                def group = meta.group
                def pipeline = meta.pipeline
                tuple(group, [(pipeline): csv])
            }.groupTuple()
            .map { group, data_list ->
                tuple(group, data_list.collectEntries())
            }

        mode_ch = channel.of("sample", "patient")

        mode_ch
            .combine(csvs_with_date)
            .map { mode, group, csv_map ->
                tuple([group: group, mode: mode], csv_map)
            }
            .set { ch_formatted_input }


        clinical_data = FORMAT_CLINICAL(
            ch_formatted_input,
            id_linking_file
        )

        meta_text = Channel.of("""cancer_study_identifier: add_text
genetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_clinical_sample.txt
        """,
        """cancer_study_identifier: add_text
genetic_alteration_type: CLINICAL
datatype: PATIENT_ATTRIBUTES
data_filename: data_clinical_patient.txt
        """)

        file_names = Channel.of("clinical_sample", "clinical_patient")

        GENERATE_META_FILE(
            file_names,
            meta_text
        )

        emit:
            csvs_with_date
}
