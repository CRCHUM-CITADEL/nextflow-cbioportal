include { FORMAT_CLINICAL } from '../../../modules/local/format_clinical'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'

workflow CLINICAL_AGGREGATE {
    take:
        filelist          // channel: [meta(group, pipeline, extraction_date), csv_path]
        genomic_subjects  // val: path to genomic subjects TSV, or "" to skip filtering

    main:

        csvs = filelist
            .map { meta, csv ->
                tuple(meta.group, [(meta.pipeline): csv])
            }
            .groupTuple()
            .map { group, data_list ->
                tuple(group, data_list.collectEntries())
            }
            .combine(genomic_subjects)
            .map { group, csv_map, gs ->
                if (gs) csv_map.genomic_subjects = gs
                tuple(group, csv_map)
            }

        all_groups = csvs.map { group, csv_map -> group }.unique()

        mode_ch = channel.of("sample", "patient")

        mode_ch
            .combine(csvs)
            .map { mode, group, csv_map ->
                return tuple([group: group, mode: mode], csv_map)
            }
            .set { ch_formatted_input }

        FORMAT_CLINICAL(ch_formatted_input)

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

        all_groups_times_two = all_groups.combine(file_names).map { g, f -> g }

        GENERATE_META_FILE(
            all_groups_times_two,
            file_names,
            meta_text
        )

    emit:
        csvs
}
