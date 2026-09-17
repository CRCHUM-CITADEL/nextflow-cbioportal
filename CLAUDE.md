# CLAUDE.md

Monorepo of three independent Nextflow DSL2 pipelines that transform sequencing/clinical data into [cBioPortal](https://www.cbioportal.org/) format. All target HPC (SLURM + Apptainer — no Docker, no sudo).

| Pipeline               | Input                                                | Purpose                                     |
| ---------------------- | ---------------------------------------------------- | ------------------------------------------- |
| `ampliseq/`            | Ampliseq VCFs + TSV exports                          | Mutations, CNAs, SVs, clinical → cBioPortal |
| `oncoanalyser_3.0.0/`  | Oncoanalyser 3.0.0 output (WGS/WTS) + clinical CSVs  | Genomic + clinical + ML → cBioPortal        |
| `dragen_4.4/`          | DRAGEN 4.4 somatic/germline + clinical CSVs          | Sister pipeline to oncoanalyser             |

A pipeline directory is named after the **upstream tool version it consumes**, not after its
own release version — `oncoanalyser_3.0.0/` reads oncoanalyser 3.0.0 output while itself
being at pipeline release 4.0.0.

Each pipeline is self-contained with its own `main.nf`, `nextflow.config`, `modules/`, `subworkflows/`, `bin/`, and `tests/`. Consult pipeline-specific `CLAUDE.md` when working inside a subdirectory.

## Shared Conventions

- nf-core DSL2 structure: `workflows/` → `subworkflows/local/` → `modules/local/`
- `bin/` scripts are R (optparse, data.table) or Python
- Containers: ORAS-hosted at `ghcr.io/crchum-citadel/`; `apptainer.cacheDir` = `containers/`
- Always add `stub:` blocks to new processes
- Use `params.container_*` — never hardcode image paths
- Output missing values as `NA`
- Site-specific settings do not belong in `nextflow.config`. It keeps portable
  defaults only (public container images, empty reference paths); CRCHUM values
  live in a sibling `nextflow_<site>.config` loaded by a matching profile. The
  oncoanalyser pipeline establishes this with `nextflow_citadel.config` and
  `-profile citadel` — list the site profile LAST so it wins on overlap.
  `nextflow.config` itself cannot be renamed: Nextflow auto-loads it, and it
  holds the `test` profile, the nf-schema plugin and the manifest.
- Nextflow optional outputs: `optional: true` is an option on the whole output
  declaration, not an argument to `path()` — `path("f.txt", optional: true)`
  silently does nothing

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
