include { FORMAT_ML_CNV } from '../../../modules/local/format_ml_cnv'
include { FORMAT_ML_EXPRESSION } from '../../../modules/local/format_ml_expression'
include { FORMAT_ML_MUTATION } from '../../../modules/local/format_ml_mutation'
<<<<<<< HEAD
include { FORMAT_PROCESS_ML_SV } from '../../../modules/local/format_process_ml_sv'
include { PROCESS_ML_CNV } from '../../../modules/local/process_ml_cnv'
include { PROCESS_ML_EXPRESSION } from '../../../modules/local/process_ml_expression'
include { PROCESS_ML_MUTATION } from '../../../modules/local/process_ml_mutation'
include { DOWNLOAD_KNOWN_FUSIONS } from '../../../modules/local/download_known_fusions'
include { DOWNLOAD_CANCER_HOTSPOTS } from '../../../modules/local/download_cancer_hotspots'

workflow GENOMIC_ML {
    take:
        cnv_result_long     // tuple (group, filepath)
        expression_result   // tuple (group, filepath)
        mutation_result     // tuple (group, filepath)
        sv_result           // tuple (group, filepath)
        cosmic_data         // file

    main:
        
        known_fusions_data = DOWNLOAD_KNOWN_FUSIONS()
        cancer_hotspot_data = DOWNLOAD_CANCER_HOTSPOTS()

        FORMAT_ML_CNV(
            cnv_result_long
        )

        PROCESS_ML_CNV(
           FORMAT_ML_CNV.out    
        )

        FORMAT_ML_EXPRESSION(
            expression_result
        )

        PROCESS_ML_EXPRESSION(
           FORMAT_ML_EXPRESSION.out
        )

        FORMAT_ML_MUTATION(
            mutation_result
        )

        PROCESS_ML_MUTATION(
            FORMAT_ML_MUTATION.out
        )

        FORMAT_PROCESS_ML_SV(
            cosmic_data,
            known_fusions,
            sv_results
        )

    emit:
        ml_sv           = FORMAT_PROCESS_ML_SV.out
    //    ml_mutation     = PROCESS_ML_MUTATION.out
    //    ml_expression   = PROCESS_ML_EXPRESSION.out
        ml_cnv          = PROCESS_ML_CNV.out
=======

workflow GENOMIC_ML {
    take:
        groups
        cnv_result_long
        expression_result
        mutation_result
        //TODO : sv results

    main:

        FORMAT_ML_CNV(
            groups,
            cnv_result_long
        )

        FORMAT_ML_EXPRESSION(
            groups,
            expression_result
        )

        FORMAT_ML_MUTATION(
            groups,
            mutation_result
        )
>>>>>>> 555063098b4c3191c53c1bb5f99da149294529b1
}
