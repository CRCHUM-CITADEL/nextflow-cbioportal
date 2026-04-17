include { ESVEE_SV_TO_CBIOPORTAL } from '../../../modules/local/esvee_sv_to_cbioportal'

workflow GENOMIC_SV {
    take:
        esvee_vcf           // tuple (meta, *-T.esvee.unfiltered.vcf.gz) — tumor only
        ensembl_annotations // path — used by ESVEE for gene-coordinate annotation

    main:
        esvee_sv = ESVEE_SV_TO_CBIOPORTAL(esvee_vcf, ensembl_annotations)

    emit:
        sv_out = esvee_sv.sv // channel [ meta, data_sv.txt ]
}
