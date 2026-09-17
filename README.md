# CRCHUM-CITADEL/nextflow-cbioportal

This monorepo contains three independent Nextflow DSL2 pipelines that transform sequencing and clinical data into [cBioPortal](https://www.cbioportal.org/)-compatible format. All pipelines target HPC environments (SLURM + Apptainer — no Docker, no sudo).

## Pipelines

| Pipeline              | Input                                               | Purpose                                           | README                                                 |
| --------------------- | --------------------------------------------------- | ------------------------------------------------- | ------------------------------------------------------ |
| `ampliseq/`           | Ampliseq VCFs + TSV exports                         | Mutations, CNAs, SVs, clinical files → cBioPortal | [ampliseq/README.md](ampliseq/README.md)               |
| `oncoanalyser_3.0.0/` | Oncoanalyser 3.0.0 output (WGS/WTS) + clinical CSVs | Genomic + clinical files + ML tables → cBioPortal | [oncoanalyser_3.0.0/README](oncoanalyser_3.0.0/README) |
| `dragen_4.4/`         | DRAGEN 4.4 somatic/germline + clinical CSVs         | Genomic + clinical files + ML tables → cBioPortal | [dragen_4.4/README.md](dragen_4.4/README.md)           |

Each directory is named after the upstream tool version it consumes, not after its own
release version.

Each pipeline is self-contained. Consult the pipeline-specific README for setup, configuration, samplesheet format, and run instructions.
