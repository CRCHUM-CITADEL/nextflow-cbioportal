//
// Subworkflow with functionality specific to the CRCHUM-CITADEL/nextflow-sante-precision pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { validateParameters; samplesheetToList } from 'plugin/nf-schema'

include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args // array: List of positional nextflow CLI args
    mode              // string: pipeline mode [clinical, genomic]
    outdir            // string: The output directory where the results will be saved
    genomic_input     // string: Path to input samplesheet
    clinical_input    // string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //

    validateInputParameters()
    validateParameters()

    //
    // Create channel from input file provided through params.input
    //
    if (mode == 'clinical'){

        samplesheet_list = Channel.fromList(samplesheetToList(clinical_input, "assets/schema_clinical_input.json"))

    } else if (mode == 'genomic'){
        if (!params.gencode_annotations){
            error("ERROR: Missing gencode_annotations file (tsv format) Check input in nextflow.config")
        }

        if (!params.ensembl_annotations){
            error("ERROR: Missing gencode_annotations file (tsv format) Check input in nextflow.config")
        }

        samplesheet_list = Channel.fromList(samplesheetToList(genomic_input, "assets/schema_genomic_input.json"))
    } else {
        error("ERROR: This should not be possible, the mode check should have caught this. Killing pipeline.")
    }


    emit:
    samplesheet = samplesheet_list
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:
    summary_params = [:]

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            // TODO: wait for HPC access
            // completionEmail(
            //     summary_params,
            //     email,
            //     email_on_fail,
            //     plaintext_email,
            //     outdir,
            //     monochrome_logs,
            //     []
            // )
        }

        completionSummary(monochrome_logs)
        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {

    // check modes and input
    if (!params.mode){
        error("ERROR: Pipeline mode not chosen in configuration file. Choices : 'genomic' or 'clinical'")
    }
    params.mode = params.mode.toLowerCase()
    if ( !params.mode in ['genomic','clinical'] ) {
        error("Error: Invalid pipeline mode chosen. Choices : 'genomic' or 'clinical'")
    }

    // make sure there's input
    if (params.mode == "genomic" && !params.genomic_samplesheet){
        error("ERROR: Could not find genomic samplesheet. Not running any tests. Check input in nextflow.config")
    }

    if (params.mode == "genomic") {
        if (!params.genome_reference) {
            error("ERROR: genome_reference parameter is required for genomic mode")
        }

        def genome_file = file(params.genome_reference)

        if (!genome_file.exists()) {
            error("ERROR: Genome reference file does not exist: ${params.genome_reference}")
        }

    }

    if (params.mode == "clinical" && !params.clinical_samplesheet){
        error("ERROR: Could not find clinical filesheet. Not running any tests. Check input in nextflow.config")
    }

    if (params.mode == "clinical") {
        if (!params.id_linking_file) {
            error("ERROR: id_linking_file parameter is required for genomic mode")
        }

        def genome_file = file(params.id_linking_file)

        if (!genome_file.exists()) {
            error("ERROR: Genome reference file does not exist: ${params.id_linking_file}")
        }

    }


}
