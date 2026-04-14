# CLAUDE.md — oncoanalyser

Nextflow DSL2 pipeline (nf-core style) that processes genomic and clinical data into cBioPortal-ready output. Runs on HPC (SLURM) using Apptainer containers — **no sudo access, no Docker**.

---

## Tech Stack

| Layer | Technology |
|---|---|
| Workflow orchestration | Nextflow DSL2 (≥24.10.2) |
| Container runtime | Apptainer (Singularity-compatible, rootless) |
| Data processing | R 4.5.1, Python 3.12 |
| HPC scheduler | SLURM |
| Schema validation | nf-schema 2.5.1 |
| Linting | prettier + pre-commit |
| Testing | nf-test 0.9.5 |

---

## Repository Layout

```
main.nf                     # Entry point — selects genomic or clinical workflow
nextflow.config             # Global params, profiles (apptainer, slurm, test, debug, gpu)
conf/
  base.config               # Process resource labels (process_single → process_high_memory)
  modules.config            # Per-module publishDir / ext.args overrides
workflows/
  genomic.nf                # Genomic pipeline (CNV, SV, expression, mutations, ML)
  clinical.nf               # Clinical pipeline
subworkflows/local/
  genomic_cnv/              # Copy number variation (Purple TSVs → SEG + CNA long)
  genomic_sv/               # Structural variants (ESVEE VCF → data_sv.txt)
  genomic_expression/       # Gene expression / TPM (Isofox CSV → TPM TSV)
  genomic_mutations/        # Mutation calling (PAVE VCFs → MAF via VEP + vcf2maf + PCGR)
  genomic_ml/               # ML data formatting (R scripts + COSMIC/ChimerDB fusions)
  genomic_aggregate_output/ # Final cBioPortal aggregation (collectFile + meta/case-list generation)
  clinical_aggregate/       # Clinical data aggregation
modules/local/              # Custom Nextflow process definitions
modules/nf-core/            # Pinned nf-core modules (bcftools, ensemblvep, vcf2maf)
bin/                        # R scripts executed inside containers
containers/                 # Apptainer .def files + cached .sif/.img images
assets/test_data/           # Test fixtures: precomputed files, oncoanalyser output, annotations
tests/
  clinical.nf.test          # Pipeline-level clinical test (snapshot-verified)
  subworkflows/             # Subworkflow-level tests (one per subworkflow)
  nextflow_subworkflow.config  # Shared config for subworkflow tests (overwrite = true)
```

---

## Pipeline Modes

The pipeline has two modes controlled by `params.mode`:

- **`genomic`** — processes VCFs, fusions, expression quantification → cBioPortal genomic files + ML-ready tables
- **`clinical`** — processes clinical CSVs (patient, diagnosis, treatment, specimen, etc.) → cBioPortal clinical files

Both modes share the same entry point (`main.nf`) and can be run independently or together.

---

## Containers (Apptainer — No Docker, No Sudo)

Container images are referenced as ORAS URLs and cached locally in `containers/`:

```
container_r       = "oras://ghcr.io/crchum-citadel/sdp-r:4.5.1"
container_python  = "oras://ghcr.io/crchum-citadel/sdp-python:3.12"
container_pcgr    = "oras://ghcr.io/sigven/pcgr:2.1.2.singularity"
container_vcf2maf = "oras://ghcr.io/crchum-citadel/vcf2maf_ensembl-vep:v1.6.22_2"
```

Nextflow + Apptainer caches ORAS images as `.img` files (not `.sif`) in `apptainer.cacheDir = "${projectDir}/containers"`. If you build a container manually with `apptainer build`, also create a `.img` symlink so Nextflow's cache lookup succeeds:

```bash
apptainer build --fakeroot containers/ghcr.io-crchum-citadel-sdp-r-4.5.1.sif containers/r_v4.5.1.def
ln -sf ghcr.io-crchum-citadel-sdp-r-4.5.1.sif containers/ghcr.io-crchum-citadel-sdp-r-4.5.1.img
```

**Container rules:**
- Always use the `apptainer` profile when running locally or on HPC: `-profile apptainer`
- Never suggest Docker commands or Docker-specific syntax
- R scripts in `bin/` run inside `container_r`; Python scripts run inside `container_python`
- `PYTHONNOUSERSITE=1`, `R_PROFILE_USER=/.Rprofile`, `R_ENVIRON_USER=/.Renviron` are set globally to prevent host library conflicts

---

## Process Resource Labels

Always use these labels in new process definitions (defined in `conf/base.config`):

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
| `process_gpu` | — | — | GPU-aware |

All values scale with `task.attempt` on retry.

---

## Key Commands

### Run tests
```bash
# All locally-runnable tests
./nf-test test tests/clinical.nf.test tests/subworkflows/ --profile test,apptainer

# Individual subworkflow test
./nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer

# Update snapshots after intentional output change
./nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer --update-snapshot
```

### Run pipeline locally (apptainer, no SLURM)
```bash
nextflow run main.nf -profile apptainer --mode genomic --genomic_samplesheet samplesheet.csv
```

### Run on HPC (SLURM + Apptainer)
```bash
nextflow run main.nf -profile slurm,apptainer --mode genomic --genomic_samplesheet samplesheet.csv
```

### Lint (run before committing)
```bash
pre-commit run --all-files
```

### Install nf-core module
```bash
nf-core modules install <module_name>
```

---

## HPC Constraints (No Sudo)

- **No `sudo`, no Docker** — Apptainer only
- **No internet on compute nodes** — set `NXF_OFFLINE=true` and pre-pull containers on login nodes
- **SLURM account** is `def-chasse` (set in `nextflow.config` `slurm` profile)
- Avoid writing to system paths; use `$SCRATCH`, `$PROJECT`, or relative paths from `projectDir`
- Large reference files (VEP cache, PCGR data) must be pre-staged; pipeline will attempt download if paths are empty — this will fail on compute nodes
- SLURM profile limits: `queueSize = 250`, `submitRateLimit = '10 sec'`, `pollInterval = '30 sec'`

---

## Writing Nextflow Modules

Follow nf-core DSL2 conventions for all new modules:

```groovy
process MY_PROCESS {
    tag "$meta.id"
    label 'process_medium'

    container params.container_r    // always use a params.container_* variable

    input:
    tuple val(meta), path(input_file)

    output:
    tuple val(meta), path("*.tsv"), emit: results

    script:
    """
    Rscript ${projectDir}/bin/my_script.R \\
        --input $input_file \\
        --output ${meta.id}.tsv
    """

    stub:
    """
    touch ${meta.id}.tsv
    """
}
```

- Place new local modules under `modules/local/<module_name>/main.nf`
- Place new local subworkflows under `subworkflows/local/<subworkflow_name>/main.nf`
- nf-core modules go under `modules/nf-core/` — never edit these manually; use `nf-core modules update`
- Use `ext.args` in `conf/modules.config` for per-module argument overrides, not hardcoded values
- Always add a `stub:` block — it allows `-stub-run` in tests to bypass heavy containers

---

## R Scripts in `bin/`

- All R scripts use `optparse` for CLI argument parsing
- Key libraries: `tidyverse`, `data.table`, `stringr`, `R.utils`
- Scripts are called from Nextflow process `script:` blocks via `Rscript ${projectDir}/bin/script.R`
- Do not add `renv` or install packages at runtime — all packages must be baked into the container `.def`
- To add a new R package: edit `containers/r_v4.5.1.def` and rebuild the container
- `ml_format_cnv.R` and `ml_format_expression.R` check `basename(input)` against `"data_cna_long.txt"` / `"data_expression.txt"` — input files must use those exact names
- **Output formatting**: empty/missing values must be `NA` (not `.` or empty string). Use `write.table(..., na = "NA")`. For `data_sv.txt`, rows must have Hugo gene symbols at both sites — filter out rows without gene annotation.

---

## Code Style

- Prettier enforces formatting (`.prettierrc.yml`) — run `pre-commit run --all-files` before pushing
- nf-core modules are excluded from linting (`--path modules/nf-core`)
- Nextflow files: 4-space indentation, single quotes for strings
- R scripts: tidyverse style, snake_case variable names
- Shell scripts in process blocks: use `\\` line continuation, `set -euo pipefail` is enforced globally via `process.shell`

---

## Testing

Tests use **nf-test** (binary at `./nf-test` in the pipeline directory) with the `nft-utils@0.0.3` plugin.

### Test structure

| Test file | Type | Containers needed |
|---|---|---|
| `tests/clinical.nf.test` | Pipeline-level | R only |
| `tests/subworkflows/clinical_aggregate.nf.test` | Subworkflow | R only |
| `tests/subworkflows/genomic_cnv.nf.test` | Subworkflow | R only |
| `tests/subworkflows/genomic_sv.nf.test` | Subworkflow | R only |
| `tests/subworkflows/genomic_expression.nf.test` | Subworkflow | R only |
| `tests/subworkflows/genomic_aggregate_output.nf.test` | Subworkflow | R only |
| `tests/subworkflows/genomic_ml.nf.test` | Subworkflow (fully stubbed) | None |
| `tests/subworkflows/genomic_mutations.nf.test` | Subworkflow | PCGR + vcf2maf (CI only) |

### Test data

- **`assets/test_data/oncoanalyser_output/TEST/T{1,2}/`** — per-subject tool outputs (PAVE VCFs, ESVEE VCFs, Purple TSVs, Isofox CSVs)
- **`assets/test_data/precomputed/`** — pre-processed files for aggregate/ML tests
- **`assets/test_data/annotations/`** — BioMart TSV (chr21, includes `entrez_ncbi_id`) + ChimerKB4.xlsx
- **`assets/test_data/clinical/`** — clinical CSVs for clinical pipeline tests

To regenerate test fixtures from scratch:
```bash
bash bin/create_test_data.sh
```

### nf-test conventions

- Channel file outputs are `String` in `then {}` blocks — use `path(f.toString())` not `f.name` or `f.readLines()`
- Sort snapshots by filename for stable ordering: `.sort { it.toString().split('/').last() }`
- `collectFile` with `storeDir` does NOT auto-create parent directories — tests must call `file("${params.outdir}/GROUP").mkdirs()` before workflow input
- `genomic_ml` uses `options "-stub-run"` to bypass `DOWNLOAD_KNOWN_FUSIONS` curl download; all ML format/process modules have `stub:` blocks

### Snapshots

Snapshots live in `tests/*.snap` and `tests/subworkflows/*.snap`. Update with `--update-snapshot` when output intentionally changes.

---

## Parameter Reference (Key Params)

| Param | Purpose |
|---|---|
| `mode` | `"genomic"` or `"clinical"` |
| `genomic_samplesheet` | CSV with columns `group, subject_id, sample_id, folder` (one row per subject) |
| `clinical_samplesheet` | CSV with clinical data paths |
| `id_linking_file` | Links genomic ↔ clinical sample IDs |
| `ensembl_annotations` | BioMart TSV for gene annotation (must include `entrez_ncbi_id` column) |
| `vep_data` | Path to pre-staged VEP cache (empty → pipeline downloads test cache) |
| `genome_reference` | GRCh38 FASTA |
| `container_r/python/vcf2maf/pcgr` | ORAS image URIs |
| `outdir` | Output directory (default: `output/`) |

---

## Branch Strategy

- `main` — stable releases
- `dev` — active development, target for PRs
- Feature branches merge into `dev` via PR; `dev` merges into `main` for releases
