# CLAUDE.md

Monorepo of three independent Nextflow DSL2 pipelines → cBioPortal, for HPC (SLURM + Apptainer; no
Docker, no sudo). Each is self-contained (`main.nf`, `nextflow.config`, `modules/`, `subworkflows/`,
`bin/`, `tests/`) and has its own `CLAUDE.md`.

| Pipeline              | Input                                            |
| --------------------- | ------------------------------------------------ |
| `ampliseq/`           | Ampliseq VCFs + TSV exports                      |
| `oncoanalyser_3.0.0/` | Oncoanalyser 3.0.0 WGS/WTS + clinical CSVs (+ML) |
| `dragen_4.4/`         | DRAGEN 4.4 somatic/germline + clinical CSVs      |

Directories are named after the upstream tool version they consume, not their own release
(`oncoanalyser_3.0.0/` is pipeline release 4.0.0).

## Shared conventions

- nf-core layout: `workflows/` → `subworkflows/local/` → `modules/local/`; `bin/` is R (optparse,
  data.table) or Python
- Containers come from `params.container_*` (never hardcoded); ORAS images at `ghcr.io/crchum-citadel/`;
  `apptainer.cacheDir = containers/`
- Every new process gets a `stub:` block. Missing values are written as `NA`
- `nextflow.config` keeps portable defaults only; site values go in `nextflow_<site>.config` behind a
  matching profile listed last (e.g. `-profile slurm,apptainer,citadel`)
- `optional: true` belongs on the output declaration (`tuple ..., emit: y, optional: true`);
  `path("f.txt", optional: true)` silently does nothing

## Commands (from the pipeline dir)

```bash
nf-test test tests/<test>.nf.test --profile test,apptainer [--update-snapshot]
pre-commit run --all-files
```

Branches: feature → `dev` via PR; `dev` → `main` for releases.
