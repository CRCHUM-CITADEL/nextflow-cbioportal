#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CRCHUM-CITADEL/nextflow-sante-precision
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/CRCHUM-CITADEL/nextflow-sante-precision
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENOMIC   } from './workflows/genomic.nf'
include { CLINICAL  } from './workflows/clinical.nf'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks and checks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.monochrome_logs,
        args,
        params.mode,
        params.outdir,
        params.genomic_samplesheet,
        params.clinical_samplesheet,
    )

    //
    // WORKFLOW: Run main workflow
    //
    // NFCORE_CITADEL_TEST (
    //     PIPELINE_INITIALISATION.out.samplesheet
    // )
    if (params.mode == 'genomic'){

        // Check params and create channels + flags at the top
        ch_vep_data = params.vep_data ? Channel.fromPath(params.vep_data) : Channel.empty()
        ch_pcgr_data = params.pcgr_data ? Channel.fromPath(params.pcgr_data) : Channel.empty()

        needs_vep_download = !params.vep_data
        needs_pcgr_download = !params.pcgr_data

        GENOMIC (
            PIPELINE_INITIALISATION.out.samplesheet,
            params.ensembl_annotations,
            params.ensembl_annotations_expr,
            ch_vep_data,
            ch_pcgr_data,
            needs_vep_download,
            needs_pcgr_download,
            params.genome_reference
        )
    }
    else if (params.mode == 'clinical'){
        CLINICAL(
            PIPELINE_INITIALISATION.out.samplesheet,
            params.id_linking_file
        )
    }


    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
