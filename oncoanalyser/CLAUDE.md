# CLAUDE.md — oncoanalyser

Nextflow DSL2 pipeline converting oncoanalyser WGS/WTS output + clinical CSVs → cBioPortal-ready files. HPC-only (SLURM + Apptainer — **no Docker, no sudo**).

---

## Layout

```
workflows/genomic.nf          # CNV, SV (DNA+RNA fusions), expression, mutations, sigs, ML
workflows/clinical.nf         # Clinical CSV processing
subworkflows/local/
  genomic_cnv/                # Purple TSVs → SEG + CNA long
  genomic_sv/                 # ESVEE tumor VCF → data_sv.txt (DNA SVs only)
  genomic_expression/         # Isofox CSV → TPM TSV
  genomic_mutations/          # PAVE VCFs → MAF via VEP + vcf2maf + PCGR
  genomic_ml/                 # ML tables (COSMIC/ChimerDB fusions)
  genomic_aggregate_output/   # collectFile merge → group-level cBioPortal files
  clinical_aggregate/
modules/local/                # One process per file
bin/                          # R scripts (optparse, data.table, tidyverse)
tests/subworkflows/           # nextflow_workflow tests
tests/modules/                # nextflow_process tests
```

---

## Containers

```
container_r            = "oras://ghcr.io/crchum-citadel/sdp-r:4.5.1"
container_python       = "oras://ghcr.io/crchum-citadel/sdp-python:3.12"
container_pcgr         = "oras://ghcr.io/sigven/pcgr:2.1.2.singularity"
container_vcf2maf      = "oras://ghcr.io/crchum-citadel/vcf2maf_ensembl-vep:v1.6.22_2"
container_sigprofiler  = "oras://ghcr.io/crchum-citadel/sdp-sigprofiler:1.1.3"
```
- Always `-profile apptainer`. Never Docker.
- R scripts → `container_r`; Python scripts → `container_python`; SigProfiler scripts → `container_sigprofiler`.
- `apptainer build` creates `.sif`; Nextflow expects `.img` — symlink after building.

---

## Process Labels (`conf/base.config`)

| Label | CPUs | Memory | Time |
|---|---|---|---|
| `process_single` | 1 | 6 GB | 4 h |
| `process_low` | 2 | 12 GB | 4 h |
| `process_medium` | 6 | 36 GB | 8 h |
| `process_high` | 12 | 72 GB | 16 h |
| `process_high_memory` | — | 200 GB | — |
| `error_retry` | — | — | maxRetries 3 |

---

## HPC Constraints

- No internet on compute nodes — set `NXF_OFFLINE=true`; pre-pull containers on login nodes.
- SLURM account: `def-chasse`.
- VEP/PCGR data must be pre-staged.

---

## Modules & R Scripts — Key Rules

- Always add a `stub:` block to every new process.
- Use `params.container_*` — never hardcode image paths.
- Output missing values as `NA`. Use `write.table(..., na = "NA")`.
- `data_sv.txt` rows require Hugo symbols at both sites — filter unannotated rows.
- **SV classification** (`gen_esvee_sv_to_cbioportal.R`): BND ALT strand → `(+,-)` DELETION, `(-,+)` DUPLICATION, `(+,+)/(−,−)` INVERSION, diff chr TRANSLOCATION. DNA SVs: `DNA_Support=Yes, RNA_Support=No`.
- **RNA fusions** (`gen_isofox_fusion_to_cbioportal.R`): reads `*.isf.pass_fusions.csv`; auto-detects columns across Isofox versions. `Class=FUSION, DNA_Support=No, RNA_Support=Yes`. Both merge into `data_sv.txt`.
- `ml_format_cnv.R` / `ml_format_expression.R` check `basename(input)` — inputs must be named `data_cna_long.txt` / `data_expression.txt`.

---

## Mutational Signatures

Three signature types are fitted against COSMIC v3.6 reference (`assets/COSMIC_Human_v3.6.zip`), producing 6 cBioPortal GENERIC_ASSAY data files per group. Metadata files in `assets/` provide `category`, `etiology`, `main_effect` for NAME/DESCRIPTION columns.

### SBS — Single Base Substitutions (101 signatures)

**Container:** `container_sigprofiler` (SigProfilerAssignment + pandas + scipy + pysam)

**Contribution** (`data_mutational_signatures_contribution_SBS.txt`):
- Input: `{subject}-T.sig.snv_counts.csv` (96-channel trinucleotide counts)
- Script: `bin/run_sigprofiler_sbs.py` → `SigProfilerAssignment.Analyzer.cosmic_fit(context_type="96")`
- Metadata: `assets/cosmic_sbs_metadata.tsv`
- Columns: `ENTITY_STABLE_ID` (`mutational_signatures_contribution_{SBS_ID}`), `NAME`, `DESCRIPTION`, `{sample}`
- Merge: `gen_merge_sigs_to_cbioportal.R` (R, outer join by ENTITY_STABLE_ID)

**Counts** (`data_mutational_signatures_counts_SBS.txt`):
- Input: same `sig.snv_counts.csv`
- Script: `bin/gen_sigs_counts_to_cbioportal.R` (context transform: `C>A_ACA` → `A[C>A]A`)
- Columns: `ENTITY_STABLE_ID` (`mutational_signatures_matrix_{5'}_{from}-{to}_{3'}`), `NAME`, `CATEGORY`, `{sample}`
- Merge: `gen_merge_sigs_counts_to_cbioportal.R`

### DBS — Doublet Base Substitutions (22 signatures)

**Contribution** (`data_mutational_signatures_contribution_DBS.txt`):
- Input: `{subject}-T.pave.somatic.vcf.gz` (extracts adjacent SNV pairs + 2bp MNVs → 78-channel DBS matrix)
- Script: `bin/run_sigprofiler_dbs.py` → `Analyzer.cosmic_fit(context_type="DINUC")`
- Metadata: `assets/cosmic_dbs_metadata.tsv`
- Strand normalization to 10 canonical ref dinucleotides (AC, AT, CC, CG, CT, GC, TA, TC, TG, TT)

**Counts** (`data_mutational_signatures_counts_DBS.txt`):
- Same script produces both contribution and counts
- `ENTITY_STABLE_ID` = `mutational_signatures_matrix_DBS_{ref}_{alt}`, `CATEGORY` = reference dinucleotide

### SV — Structural Variants (12 signatures)

**Contribution** (`data_mutational_signatures_contribution_SV.txt`):
- Input: `{subject}-T.esvee.somatic.vcf.gz` (BND records → 32-channel SV type matrix)
- Script: `bin/run_sigprofiler_sv.py` → `scipy.optimize.nnls` (SigProfilerAssignment doesn't support SV context)
- Metadata: `assets/cosmic_sv_metadata.tsv`
- SV classification: type (del/inv/tds/trans) × size bin (5 bins) × clustering (25kb proximity threshold) = 32 categories

**Counts** (`data_mutational_signatures_counts_SV.txt`):
- Same script produces both contribution and counts
- `ENTITY_STABLE_ID` = `mutational_signatures_matrix_SV_{type}`, `CATEGORY` = base type (DEL/INV/TDS/TRANS)

### Common

- All meta files use `generic_assay_type: MUTATIONAL_SIGNATURE`, `datatype: LIMIT-VALUE`, `pivot_threshold_value: 0.0`.
- Merge modules (6 total) reuse `gen_merge_sigs_to_cbioportal.R` / `gen_merge_sigs_counts_to_cbioportal.R` with different output filenames.
- Zero-mutation edge case: scripts write header-only files; <50 mutations: stderr warning but continues.

---

## Testing

```bash
# All locally-runnable tests
./nf-test test tests/clinical.nf.test tests/subworkflows/ tests/modules/ --profile test,apptainer

# Single test + update snapshot
./nf-test test tests/subworkflows/genomic_sv.nf.test --profile test,apptainer --update-snapshot
```

| Test file | Containers |
|---|---|
| `tests/subworkflows/genomic_cnv.nf.test` | R only |
| `tests/subworkflows/genomic_sv.nf.test` | R only |
| `tests/subworkflows/genomic_expression.nf.test` | R only |
| `tests/subworkflows/genomic_aggregate_output.nf.test` | R only |
| `tests/subworkflows/genomic_ml.nf.test` | None (fully stubbed) |
| `tests/subworkflows/genomic_mutations.nf.test` | PCGR + vcf2maf (CI only) |
| `tests/modules/isofox_fusion_to_cbioportal.nf.test` | R only |
| `tests/modules/sigprofiler_sbs.nf.test` | SigProfiler |
| `tests/modules/sigprofiler_dbs.nf.test` | SigProfiler |
| `tests/modules/sigprofiler_sv.nf.test` | SigProfiler |
| `tests/modules/sigs_counts_to_cbioportal.nf.test` | R only |

**nf-test gotchas:**
- Use `path(f.toString())` — channel file outputs are `String`, not `Path`.
- Sort snapshots: `.sort { it.toString().split('/').last() }`.
- `collectFile` with `storeDir` won't create directories — call `file("${params.outdir}/GROUP").mkdirs()` in test setup.
- `genomic_ml` uses `options "-stub-run"` to skip the `DOWNLOAD_KNOWN_FUSIONS` curl call.
- Module tests use `nextflow_process {}` blocks; snapshots live in `tests/modules/*.snap`.
