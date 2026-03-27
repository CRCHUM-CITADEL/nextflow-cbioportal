include { VCF2MAF           } from '../../../modules/local/vcf2maf'
include { DOWNLOAD_VEP_TEST } from '../../../modules/local/download_vep_test'
include { GENERATE_CASE_LIST} from '../../../modules/local/generate_case_list'
include { GENERATE_META_FILE} from '../../../modules/local/generate_meta_file'

workflow GENOMIC_MUTATIONS {
    take:
        som_dna_vcf  // tuple (meta, sage.vcf.gz) — SAGE somatic tumor-normal VCF
        germ_dna_vcf // tuple (meta, sage.germline.vcf.gz) — SAGE germline VCF (normal sample)
        fasta        // path — GRCh38 reference FASTA
        vep_data     // channel<path> — pre-staged VEP cache directory
        needs_vep    // boolean — download VEP test cache if true

    main:
        ch_vep_data = needs_vep ? DOWNLOAD_VEP_TEST().cache_dir.first() : vep_data.first()

        VCF2MAF(som_dna_vcf.mix(germ_dna_vcf), fasta, ch_vep_data)

    emit:
        out = VCF2MAF.out.maf // channel [ meta, sample.maf ]
}
