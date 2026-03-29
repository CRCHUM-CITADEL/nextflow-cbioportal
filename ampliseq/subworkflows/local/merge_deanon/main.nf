include { MERGE_SV         } from '../../../modules/local/merge_sv/main'
include { MERGE_CNA        } from '../../../modules/local/merge_cna/main'
include { MERGE_MUTATIONS  } from '../../../modules/local/merge_mutations/main'
include { DEANON_MUTATIONS } from '../../../modules/local/deanon_mutations/main'
include { DEANON_SV        } from '../../../modules/local/deanon_sv/main'
include { DEANON_CNA       } from '../../../modules/local/deanon_cna/main'

workflow MERGE_DEANON {

    take:
    ch_sv_files       // channel: collected per-sample sv files
    ch_cna_files      // channel: collected per-sample cna files
    ch_mutation_files // channel: collected per-sample mutation files
    ch_linking        // value channel: linking file

    main:
    MERGE_SV(ch_sv_files)
    MERGE_CNA(ch_cna_files)
    MERGE_MUTATIONS(ch_mutation_files)

    DEANON_MUTATIONS(MERGE_MUTATIONS.out, ch_linking)
    DEANON_SV(MERGE_SV.out, ch_linking)
    DEANON_CNA(MERGE_CNA.out, ch_linking)

    emit:
    mutations = DEANON_MUTATIONS.out
    sv        = DEANON_SV.out
    cna       = DEANON_CNA.out
}
