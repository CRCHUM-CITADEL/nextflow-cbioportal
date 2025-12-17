include { FORMAT_ML_CNV } from '../../../modules/local/format_ml_cnv'
include { FORMAT_ML_EXPRESSION } from '../../../modules/local/format_ml_expression'
include { FORMAT_ML_MUTATION } from '../../../modules/local/format_ml_mutation'

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
}