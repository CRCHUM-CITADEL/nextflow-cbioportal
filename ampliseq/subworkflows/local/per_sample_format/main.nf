include { FORMAT_SV             } from '../../../modules/local/format_sv/main.nf'
include { FORMAT_CNA            } from '../../../modules/local/format_cna/main.nf'
include { STUB_MAF              } from '../../../modules/local/stub_maf/main.nf'
include { VCF_TO_MAF            } from '../../../modules/local/vcf_to_maf/main.nf'
include { FILTER_MUTATIONS      } from '../../../modules/local/filter_mutations/main.nf'
include { PASSTHROUGH_MUTATIONS } from '../../../modules/local/passthrough_mutations/main.nf'
include { VCF_TO_SEG            } from '../../../modules/local/vcf_to_seg/main.nf'

workflow PER_SAMPLE_FORMAT {

    take:
    ch_tsv       // channel: tuple(meta, tsv)
    ch_vcf_input // channel: tuple(meta, sample_folder)

    main:
    ch_sv_input = ch_tsv.join(ch_vcf_input)  // → tuple(meta, tsv, sample_folder)
    FORMAT_SV(ch_sv_input)
    FORMAT_CNA(ch_tsv)
    ch_sv = FORMAT_SV.out
    ch_cna = FORMAT_CNA.out

    // -------------------------------------------------------------------------
    // Segmentation: CNV VCF → .seg
    // -------------------------------------------------------------------------
    VCF_TO_SEG(ch_vcf_input)
    ch_seg = VCF_TO_SEG.out

    // -------------------------------------------------------------------------
    // Mutations: VCF → MAF, then filter by TSV coordinates
    // -------------------------------------------------------------------------
    if (params.skip_vcf2maf) {
        STUB_MAF(ch_vcf_input)
        ch_maf = STUB_MAF.out
    } else {
        if (!params.vcf2maf_container) {
            error "params.vcf2maf_container must be set when skip_vcf2maf is false"
        }
        if (!params.ref_fasta) {
            error "ref_fasta must be set when skip_vcf2maf is false"
        }
        if (!params.vep_data) {
            error "params.vep_data must be set when skip_vcf2maf is false"
        }
        ch_vep       = Channel.value(file(params.vep_data))
        ch_ref_fasta = Channel.value(file(params.ref_fasta))
        VCF_TO_MAF(ch_vcf_input, ch_vep, ch_ref_fasta)
        ch_maf = VCF_TO_MAF.out
    }

    if (params.filter_tsv_variants) {
        FILTER_MUTATIONS(ch_maf.join(ch_tsv))
        ch_mutations = FILTER_MUTATIONS.out
    } else {
        PASSTHROUGH_MUTATIONS(ch_maf)
        ch_mutations = PASSTHROUGH_MUTATIONS.out
    }

    emit:
    sv        = ch_sv
    cna       = ch_cna
    mutations = ch_mutations
    seg       = ch_seg
}
