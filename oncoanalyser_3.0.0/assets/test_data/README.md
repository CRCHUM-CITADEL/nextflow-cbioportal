# Test data (oncoanalyser_3.0.0)

Minimal fixtures for nf-test and the Python tests. Everything is chr21-scale and synthetic.

| Path                                  | What it is                                                                                                                                                                                                                                                     | Made by                                                               |
| ------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------- |
| `oncoanalyser_output/TEST/T{1,2}/`    | Oncoanalyser 3.0 output per subject: `purple/` (somatic + germline + SV VCFs, CNV, purity), `chord/` (T1 only, so T2 exercises HRD = `NA`), `isofox/`, `sage_append/`, `sigs/`. T1's SV VCF also has PON, INFERRED, single-breakend and TINC_RECOVERED records | `bin/create_test_data.sh` (`sigs/` is hand-made)                      |
| `oncoanalyser_output/TEST/T3/purple/` | Duplicate-coordinate BND pairs for the SV dedupe test                                                                                                                                                                                                          | `bin/create_test_data.sh`                                             |
| `precomputed/`                        | Per-subject outputs a real run publishes, staged as the incremental cache. `T1`/`T2` names feed subworkflow tests; `T1-1DT`/`T2-T` names feed pipeline tests                                                                                                   | `create_test_data.sh` for `T1`/`T2` + group files; the rest hand-made |
| `annotations/`                        | BioMart gene table (chr21 + test genes), COSMIC fusion sample                                                                                                                                                                                                  | hand-made                                                             |
| `reference/`                          | chr21 FASTA for SigProfiler ID                                                                                                                                                                                                                                 | downloaded                                                            |
| `resources/`                          | Mini resource kit used by the `test` profile (`resources_dir`)                                                                                                                                                                                                 | hand-made                                                             |
| `clinical/`, `samplesheet/`           | MOHCCN-style clinical CSVs and the genomic/clinical samplesheets                                                                                                                                                                                               | hand-made                                                             |

Regenerate with `bash bin/create_test_data.sh` (needs `bgzip` and `tabix`). It reproduces the
committed files exactly, so `git status` should show no change; if it does, the script and the
fixtures have drifted.

The samplesheets' `folder` paths are relative to the pipeline directory. Under nf-test the launch
directory is elsewhere, so pipeline tests that need the oncoanalyser folders write an absolute
samplesheet outside `outputDir` in `setup` (see `tests/both.nf.test`, `tests/incremental_pipeline.nf.test`).
