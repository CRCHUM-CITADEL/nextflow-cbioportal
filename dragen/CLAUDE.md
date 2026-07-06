# CLAUDE.md — dragen

Nextflow DSL2 pipeline converting DRAGEN somatic/germline + RNA output + clinical CSVs → cBioPortal-ready files. Sister pipeline to `oncoanalyser/` — same structure and conventions. HPC-only (SLURM + Apptainer — **no Docker, no sudo**).

---

## Layout

```
workflows/genomic.nf          # CNV, SV, expression, mutations, ML
workflows/clinical.nf         # Clinical CSV processing
subworkflows/local/
  genomic_cnv/                # Copy-number → SEG + CNA long
  genomic_sv/                 # SV VCFs → data_sv.txt
  genomic_expression/         # RNA quant → TPM TSV
  genomic_mutations/          # VCFs → MAF via VEP + vcf2maf + PCGR
  genomic_ml/                 # ML tables (COSMIC/ChimerDB fusions)
  genomic_aggregate_output/   # collectFile merge → group-level cBioPortal files
  clinical_aggregate/
modules/local/                # One process per file
bin/                          # R scripts (optparse, data.table, tidyverse)
tests/subworkflows/           # nextflow_workflow tests (one per subworkflow)
tests/modules/                # nextflow_process tests
```

---

## Containers

Images are ORAS-hosted, cached as `.img` in `containers/` (`apptainer.cacheDir`):
```
container_r       = "oras://ghcr.io/crchum-citadel/sdp-r:4.5.1"
container_python  = "oras://ghcr.io/crchum-citadel/sdp-python:3.12"
container_pcgr    = "oras://ghcr.io/sigven/pcgr:2.1.2.singularity"
container_vcf2maf = "oras://ghcr.io/crchum-citadel/vcf2maf_ensembl-vep:v1.6.22_2"
```
- Always use `-profile apptainer`. Never suggest Docker.
- `apptainer build` creates `.sif`; Nextflow expects `.img` — create a symlink after building manually.
- R scripts → `container_r`; Python scripts → `container_python`.

---

## Process Labels (`conf/base.config`)

| Label | CPUs | Memory | Time |
|---|---|---|---|
| `process_single` | 1 | 6 GB | 4 h |
| `process_low` | 2 | 12 GB | 4 h |
| `process_medium` | 6 | 36 GB | 8 h |
| `process_high` | 12 | 72 GB | 16 h |
| `process_long` | — | — | 20 h |
| `process_medium_memory` | 8 | 30 GB | 2 h |
| `process_high_memory` | — | 200 GB | — |
| `error_retry` | — | — | maxRetries 3 |
| `error_ignore` | — | — | errorStrategy ignore |
| `process_gpu` | — | — | GPU accelerator when `-profile gpu` |

---

## HPC Constraints

- **No internet on compute nodes** — set `NXF_OFFLINE=true`; pre-pull containers on login nodes.
- **SLURM account**: `def-chasse` (in `nextflow.config`).
- VEP cache and PCGR data must be pre-staged; empty paths trigger download and will fail on compute nodes.

---

## Modules & R Scripts — Key Rules

- Always add a `stub:` block to every new process (`-stub-run` is used in tests).
- Use `params.container_*` variables — never hardcode image paths.
- Output empty/missing values as `NA` (not `.`). Use `write.table(..., na = "NA")`.
- `data_sv.txt` rows may have `NA` Hugo symbols — the dragen fusion script (`gen_format_dragen_fusion.R`) does not filter unannotated rows (unlike oncoanalyser which does).
- `ml_format_cnv.R` and `ml_format_expression.R` check `basename(input)` — inputs must be named `data_cna_long.txt` / `data_expression.txt`.

---

## Testing

```bash
# All locally-runnable tests
./nf-test test tests/clinical.nf.test tests/genomic.nf.test tests/subworkflows/ tests/modules/ --profile test,apptainer

# Single test + update snapshot
./nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer --update-snapshot
```

| Test file | Type |
|---|---|
| `tests/clinical.nf.test` | Workflow |
| `tests/genomic.nf.test` | Workflow |
| `tests/subworkflows/genomic_cnv.nf.test` | Subworkflow |
| `tests/subworkflows/genomic_sv.nf.test` | Subworkflow |
| `tests/subworkflows/genomic_expression.nf.test` | Subworkflow |
| `tests/subworkflows/genomic_aggregate_output.nf.test` | Subworkflow |
| `tests/subworkflows/clinical_aggregate.nf.test` | Subworkflow |
| `tests/subworkflows/genomic_ml.nf.test` | Subworkflow (stub) |
| `tests/subworkflows/genomic_mutations.nf.test` | Subworkflow (CI only) |
| `tests/modules/package_cbioportal.nf.test` | Module |

**nf-test gotchas:**
- Use `path(f.toString())` — channel file outputs are `String`, not `Path`.
- Sort snapshots: `.sort { it.toString().split('/').last() }`.
- `collectFile` with `storeDir` won't create directories — call `file("${params.outdir}/GROUP").mkdirs()` in test setup.
- `genomic_ml` uses `options "-stub-run"` to skip the `DOWNLOAD_KNOWN_FUSIONS` curl call.
- Snapshots: `tests/*.snap` and `tests/subworkflows/*.snap`.
