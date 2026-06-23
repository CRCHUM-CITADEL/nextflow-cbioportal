/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { CLINICAL_AGGREGATE } from '../subworkflows/local/clinical_aggregate'
include { GENERATE_META_FILE } from '../modules/local/generate_meta_file'
include { GENERATE_CANCER_TYPE_FILE } from '../modules/local/generate_cancer_type_file'
include { GENERATE_CLINICAL_TEMPLATE } from '../modules/local/generate_clinical_template'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLINICAL {

    take:
        file_list        // clinical samplesheet rows (may be empty channel)
        id_linking_file

    main:
        ch_versions = Channel.empty()

        if (params.clinical_samplesheet) {
            // ── Full clinical processing from clinical CSVs ──────────────────
            ch_file_list = file_list
                .map { row ->
                    def group = row[0].group
                    def sub_file = row[0].file
                    def file = sub_file.startsWith("/") ? sub_file : "${projectDir}/${sub_file}"
                    def pipeline = row[0].pipeline
                    def extraction_date = row[0].date
                    return tuple([group: group, pipeline: pipeline, extraction_date: extraction_date], file)
                }

            CLINICAL_AGGREGATE(
                ch_file_list,
                id_linking_file
            )
        } else {
            // ── Template fallback: generate minimal clinical files from the linking file ──
            log.info "No clinical samplesheet provided. Generating template clinical files from linking file."

            // Read subject_id\tsample_id lines from the linking file
            ch_linking = Channel.fromPath(id_linking_file)
                .splitCsv(header: true, sep: '\t')
                .map { row -> [subject: row.subject_id, sample: row.sample_id] }

            // Build newline-separated "subject\tsample" lines
            sample_lines = ch_linking
                .map { meta -> "${meta.subject}\t${meta.sample}" }
                .collect()
                .map { lines -> lines.join('\n') }

            // Use project name as the group
            group = params.project_name

            GENERATE_CLINICAL_TEMPLATE(
                group,
                sample_lines,
                params.icd_code ?: ""
            )

            GENERATE_CANCER_TYPE_FILE(
                channel.of(group),
                params.oncotree_code
            )

            // Generate meta files for clinical sample and cancer type
            meta_text = Channel.of(
                """cancer_study_identifier: add_text
genetic_alteration_type: CLINICAL
datatype: SAMPLE_ATTRIBUTES
data_filename: data_clinical_sample.txt
                """,
                """genetic_alteration_type: CANCER_TYPE
datatype: CANCER_TYPE
data_filename: cancer_type.txt
                """)

            file_names = Channel.of("clinical_sample", "cancer_type")

            all_groups = channel.of(group)
            all_groups_times_two = all_groups.combine(file_names).map { g, name -> g }

            GENERATE_META_FILE(
                all_groups_times_two,
                file_names,
                meta_text
            )
        }

        //
        // TASK: Aggregate software versions
        //
        softwareVersionsToYAML(ch_versions)
            .collectFile(
                storeDir: "${params.outdir}/pipeline_info",
                name: 'software_versions.yml',
                sort: true,
                newLine: true
            )

}
