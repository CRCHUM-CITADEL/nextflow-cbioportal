include { FORMAT_SV             } from '../../../modules/local/format_sv/main.nf'
include { FORMAT_CNA            } from '../../../modules/local/format_cna/main.nf'
include { STUB_MAF              } from '../../../modules/local/stub_maf/main.nf'
include { MAFSMITH              } from '../../../modules/local/mafsmith/main.nf'
include { DOWNLOAD_MAFSMITH     } from '../../../modules/local/download_mafsmith/main.nf'
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
    // Mutations: VCF → MAF, then optionally keep only mutations overlapping TSV regions
    // -------------------------------------------------------------------------
    // Nextflow 26 passes CLI values as strings (--skip_vcf2maf false → "false", which is truthy)
    if (params.skip_vcf2maf.toString().toBoolean()) {
        STUB_MAF(ch_vcf_input)
        ch_maf = STUB_MAF.out
    } else {
        if (!params.mafsmith_container) {
            error "params.mafsmith_container must be set when skip_vcf2maf is false"
        }
        // Unset, empty or a literal "null" (from --mafsmith_data null) all mean: fetch the bundle
        def mafsmith_data = params.mafsmith_data?.toString()?.trim()
        ch_mafsmith_data = mafsmith_data && mafsmith_data != 'null'
            ? Channel.value(file(mafsmith_data, checkIfExists: true))
            : DOWNLOAD_MAFSMITH(ch_vcf_input.first().map { true }).data_dir.first()
        MAFSMITH(ch_vcf_input, ch_mafsmith_data)
        ch_maf = MAFSMITH.out
    }

    if (params.filter_tsv_variants.toString().toBoolean()) {
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
