include { ESVEE_SV_TO_CBIOPORTAL      } from '../../../modules/local/esvee_sv_to_cbioportal'
include { ISOFOX_FUSION_TO_CBIOPORTAL } from '../../../modules/local/isofox_fusion_to_cbioportal'

workflow GENOMIC_SV {
    take:
        esvee_vcf           // tuple (meta, esvee.somatic.vcf.gz)  — DNA structural variants
        isofox_fusion       // tuple (meta, isofox.fusion.tsv)      — RNA fusions
        ensembl_annotations // path — used by ESVEE for gene-coordinate annotation

    main:
        esvee_sv   = ESVEE_SV_TO_CBIOPORTAL(esvee_vcf, ensembl_annotations)
        isofox_sv  = ISOFOX_FUSION_TO_CBIOPORTAL(isofox_fusion)

        // Merge DNA SVs and RNA fusions into one channel; both feed data_sv.txt
        sv_combined = esvee_sv.sv.mix(isofox_sv.sv)

    emit:
        sv_out = sv_combined // channel [ meta, data_sv.txt ]
}
