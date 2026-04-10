// include modules
include { DRAGEN_FUSION_SV_TO_CBIOPORTAL } from '../../../modules/local/dragen_fusion_sv_to_cbioportal'
include { GENERATE_META_FILE } from '../../../modules/local/generate_meta_file'

workflow GENOMIC_SV {
    take:
        sv_vcf

    main:

        cbioportal_genomic_sv_files = DRAGEN_FUSION_SV_TO_CBIOPORTAL(
            sv_vcf
        )

    emit:
        sv_out = cbioportal_genomic_sv_files
}
