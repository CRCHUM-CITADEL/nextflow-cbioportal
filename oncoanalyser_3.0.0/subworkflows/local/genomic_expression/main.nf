include { ISOFOX_EXPRESSION_TO_CBIOPORTAL } from '../../../modules/local/isofox_expression_to_cbioportal'

workflow GENOMIC_EXPRESSION {
    take:
        isofox_exp           // tuple (meta, isofox.exp.tsv)
        ensembl_annotations  // path — for Ensembl ID → Entrez ID mapping

    main:
        tpm_file_ch = ISOFOX_EXPRESSION_TO_CBIOPORTAL(isofox_exp, ensembl_annotations)

    emit:
        out = tpm_file_ch.tpm // channel [ meta, tpm.tsv ]
}
