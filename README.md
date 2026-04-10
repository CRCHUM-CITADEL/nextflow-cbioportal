# CRCHUM-CITADEL/nextflow-cbioportal

This monorepo contains three independent Nextflow DSL2 pipelines that transform sequencing and clinical data into [cBioPortal](https://www.cbioportal.org/)-compatible format. All pipelines target HPC environments (SLURM + Apptainer — no Docker, no sudo).

## Pipelines

| Pipeline | Input | Purpose | README |
| -------- | ----- | ------- | ------ |
| `ampliseq/` | Ampliseq VCFs + TSV exports | Mutations, CNAs, SVs, clinical files → cBioPortal | [ampliseq/README.md](ampliseq/README.md) |
| `oncoanalyser/` | Oncoanalyser output (WGS/WTS) + clinical CSVs | Genomic + clinical files + ML tables → cBioPortal | [oncoanalyser/README.md](oncoanalyser/README.md) |
| `dragen/` | DRAGEN somatic/germline + clinical CSVs | Genomic + clinical files + ML tables → cBioPortal | [dragen/README.md](dragen/README.md) |

Each pipeline is self-contained. Consult the pipeline-specific README for setup, configuration, samplesheet format, and run instructions.
