/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { CLINICAL_AGGREGATE } from '../subworkflows/local/clinical_aggregate'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CLINICAL {

    take:
        file_list
        id_linking_file

    main:
        ch_versions = Channel.empty()

        ch_file_list = file_list
            .map { row ->
                def group = row[0].group
                def file = "${params.input_dir}/${row[0].file}"
                def pipeline = row[0].pipeline
                def extraction_date = row[0].date
                return tuple([group: group, pipeline: pipeline, extraction_date: extraction_date], file)
            }
        

        CLINICAL_AGGREGATE(
            ch_file_list,
            id_linking_file
        )

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
