# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

This is a monorepo containing three independent Nextflow DSL2 pipelines that transform sequencing and clinical data into [cBioPortal](https://www.cbioportal.org/)-compatible format. All pipelines target HPC environments (SLURM + Apptainer — no Docker, no sudo).

| Pipeline | Input | Purpose |
|---|---|---|
| `ampliseq/` | Ampliseq VCFs + TSV exports | Mutations, CNAs, SVs, clinical files → cBioPortal |
| `oncoanalyser/` | Oncoanalyser/DRAGEN output (WGS/WTS) + clinical CSVs | Genomic + clinical files + ML tables → cBioPortal |
| `dragen/` | DRAGEN somatic/germline + clinical CSVs | Genomic + clinical files + ML tables → cBioPortal (sister pipeline to oncoanalyser) |

Each pipeline is self-contained with its own `main.nf`, `nextflow.config`, `modules/local/`, `subworkflows/local/`, `bin/`, and `tests/`. **Consult the pipeline-specific `CLAUDE.md` when working inside a subdirectory.**

---

## Shared Architecture

All pipelines follow nf-core DSL2 conventions. The canonical structure (used by oncoanalyser and dragen) is:

```
<pipeline>/
  main.nf                   # Entry point
  nextflow.config           # Params, profiles (apptainer, slurm, test, debug, gpu)
  conf/
    base.config             # Process resource labels (process_single → process_high_memory)
    modules.config          # Per-module publishDir / ext.args overrides
  workflows/                # Top-level workflow DAG
  subworkflows/local/       # Grouped process workflows
  modules/local/            # Individual process definitions
  modules/nf-core/          # Pinned nf-core modules — never edit manually
  bin/                      # R and Python scripts called from process script blocks
  containers/               # Apptainer .def files and cached .sif images
  assets/test_data/         # Test samplesheets and reference files
  tests/                    # nf-test definitions and snapshots
```

**ampliseq differences:** no `conf/` directory (config is inline in `nextflow.config`), no `modules/nf-core/`, and test data lives directly under `assets/` (no `test_data/` subdirectory).

Container images for oncoanalyser and dragen are ORAS-hosted at `ghcr.io/crchum-citadel/`. ampliseq uses a locally built `.sif` and Wave-hosted images. The `apptainer.cacheDir` points to `containers/` within each pipeline directory.

---

## Common Commands

All commands are run from within the relevant pipeline subdirectory (e.g., `cd oncoanalyser/`).

```bash
# Run tests (oncoanalyser — has clinical, clinical_template, incremental tests)
nf-test test tests/clinical.nf.test --profile test,apptainer

# Run tests (dragen — has genomic and clinical tests)
nf-test test tests/genomic.nf.test --profile test,apptainer
nf-test test tests/clinical.nf.test --profile test,apptainer

# Run tests (ampliseq — module-level tests)
nf-test test tests/modules/package_cbioportal.nf.test --profile test,apptainer

# Lint before committing
pre-commit run --all-files

# Run oncoanalyser/dragen locally (no SLURM)
nextflow run main.nf -profile apptainer --mode genomic --genomic_samplesheet samplesheet.csv

# Run oncoanalyser/dragen on HPC
nextflow run main.nf -profile slurm,apptainer --mode genomic --genomic_samplesheet samplesheet.csv

# Run ampliseq (local executor only, no --mode param)
nextflow run main.nf -profile apptainer --input samplesheet.csv

# Install / update nf-core modules (oncoanalyser/dragen only)
nf-core modules install <module_name>
nf-core subworkflows install <subworkflow_name>
```

---

## Branch Strategy

- `main` — stable releases
- `dev` — active development; target for all PRs
- Feature branches merge into `dev` via PR; `dev` merges into `main` for releases

---

## Pipeline-Specific Guidance

- **`ampliseq/CLAUDE.md`** — samplesheet format, standalone `bin/` scripts, data flow, VCF→MAF options
- **`oncoanalyser/CLAUDE.md`** — dual genomic/clinical modes, container rules, process resource labels, R script conventions, HPC constraints
- **`dragen/CLAUDE.md`** — identical structure to oncoanalyser; same conventions apply
