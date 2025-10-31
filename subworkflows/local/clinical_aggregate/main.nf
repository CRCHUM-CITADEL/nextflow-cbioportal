include { FORMAT_CLINICAL } from '../../../modules/local/format_clinical'
include { ASSIGN_DATE } from '../../../modules/local/assign_date'

workflow CLINICAL_AGGREGATE {
    take:
        filelist

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
            ch_formatted_input
        )

        

        emit:
            csvs_with_date
}
