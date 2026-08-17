# CLAUDE.md

Monorepo of three independent Nextflow DSL2 pipelines that transform sequencing/clinical data into [cBioPortal](https://www.cbioportal.org/) format. All target HPC (SLURM + Apptainer — no Docker, no sudo).

| Pipeline        | Input                                                | Purpose                                     |
| --------------- | ---------------------------------------------------- | ------------------------------------------- |
| `ampliseq/`     | Ampliseq VCFs + TSV exports                          | Mutations, CNAs, SVs, clinical → cBioPortal |
| `oncoanalyser/` | Oncoanalyser/DRAGEN output (WGS/WTS) + clinical CSVs | Genomic + clinical + ML → cBioPortal        |
| `dragen/`       | DRAGEN somatic/germline + clinical CSVs              | Sister pipeline to oncoanalyser             |

Each pipeline is self-contained with its own `main.nf`, `nextflow.config`, `modules/`, `subworkflows/`, `bin/`, and `tests/`. Consult pipeline-specific `CLAUDE.md` when working inside a subdirectory.

## Shared Conventions

- nf-core DSL2 structure: `workflows/` → `subworkflows/local/` → `modules/local/`
- `bin/` scripts are R (optparse, data.table) or Python
- Containers: ORAS-hosted at `ghcr.io/crchum-citadel/`; `apptainer.cacheDir` = `containers/`
- Always add `stub:` blocks to new processes
- Use `params.container_*` — never hardcode image paths
- Output missing values as `NA`

## Commands

```bash
# Tests (run from pipeline dir)
nf-test test tests/<test>.nf.test --profile test,apptainer
nf-test test tests/<test>.nf.test --profile test,apptainer --update-snapshot

# Lint
pre-commit run --all-files
```

## Branch Strategy

`main` = stable, `dev` = active development. Feature branches → `dev` via PR; `dev` → `main` for releases.
