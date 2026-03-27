include { PURPLE_CNV_TO_CBIOPORTAL } from '../../../modules/local/purple_cnv_to_cbioportal'

workflow GENOMIC_CNV {
    take:
        purple_cnv          // tuple (meta, purple_cnv_somatic.tsv, purple_cnv_gene.tsv)
        ensembl_annotations // path

    main:
        cbioportal_cnv_files = PURPLE_CNV_TO_CBIOPORTAL(purple_cnv, ensembl_annotations)

    emit:
        segfile  = cbioportal_cnv_files.seg  // channel [ meta, seg_file  ]
        longfile = cbioportal_cnv_files.long // channel [ meta, long_file ]
}
