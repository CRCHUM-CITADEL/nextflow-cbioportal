include { PURPLE_SV_TO_CBIOPORTAL } from '../../../modules/local/purple_sv_to_cbioportal'

workflow GENOMIC_SV {
    take:
        purple_sv_vcf       // tuple (meta, <subject>-T.purple.sv.vcf.gz) — somatic
        ensembl_annotations // path — gene overlap at each breakend

    main:
        purple_sv = PURPLE_SV_TO_CBIOPORTAL(purple_sv_vcf, ensembl_annotations)

    emit:
        sv_out = purple_sv.sv // channel [ meta, data_sv.txt ]
}
