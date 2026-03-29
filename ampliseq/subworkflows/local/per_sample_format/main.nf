include { FORMAT_SV             } from '../../../modules/local/format_sv/main'
include { FORMAT_CNA            } from '../../../modules/local/format_cna/main'
include { STUB_MAF              } from '../../../modules/local/stub_maf/main'
include { VCF_TO_MAF            } from '../../../modules/local/vcf_to_maf/main'
include { FILTER_MUTATIONS      } from '../../../modules/local/filter_mutations/main'
include { PASSTHROUGH_MUTATIONS } from '../../../modules/local/passthrough_mutations/main'

workflow PER_SAMPLE_FORMAT {

    take:
    ch_tsv       // channel: tuple(meta, tsv)
    ch_vcf_input // channel: tuple(meta, sample_folder)

    main:

    // -------------------------------------------------------------------------
    // SV and CNA: from TSV
    // -------------------------------------------------------------------------
    FORMAT_SV(ch_tsv)
    FORMAT_CNA(ch_tsv)
    ch_sv  = FORMAT_SV.out
    ch_cna = FORMAT_CNA.out

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
}
