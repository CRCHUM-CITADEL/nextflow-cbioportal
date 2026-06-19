# User's Guide

## Introduction

This pipeline converts [nf-core/oncoanalyser](https://github.com/nf-core/oncoanalyser) WGS/WTS output and clinical CSVs into [cBioPortal](https://www.cbioportal.org/)-compatible import packages. It runs on HPC clusters using SLURM and Apptainer (no Docker, no sudo).

The pipeline has two independent modes:

- **Genomic mode** — processes per-subject oncoanalyser output (mutations, copy-number, structural variants, gene expression, RNA fusions, mutational signatures) and produces merged, group-level cBioPortal data files packaged as a `.tar.gz`.
- **Clinical mode** — transforms clinical CSVs (patient demographics, diagnoses, treatments, specimens, etc.) into cBioPortal patient and sample attribute files.

### Genomic mode processing steps

| Step | Input | Output |
|---|---|---|
| Copy-number variants | PURPLE somatic + gene TSVs | `data_cna_hg38.seg`, `data_cna_long.txt` |
| Structural variants (DNA) | ESVEE somatic VCF | `data_sv.txt` (DNA rows) |
| RNA fusions | Isofox pass_fusions CSV | `data_sv.txt` (RNA rows) |
| Gene expression | Isofox gene_data CSV | `data_expression.txt` |
| Mutations | PAVE somatic + germline + RNA VCFs | `data_mutations_dna_rna_germline.txt` (MAF) |
| Mutational signatures (SBS) | SNV trinucleotide counts | Contribution + counts matrices |
| Mutational signatures (DBS) | PAVE somatic VCF | Contribution + counts matrices |
| Mutational signatures (ID) | PAVE somatic VCF | Contribution + counts matrices |
| ML feature tables | Merged group-level files | ML-formatted CNV, expression, mutation, SV |

DNA SVs and RNA fusions share the same `data_sv.txt` file; they are distinguished by `DNA_Support` / `RNA_Support` flags and the `Class` column (`DELETION`/`DUPLICATION`/`INVERSION`/`TRANSLOCATION` for DNA, `FUSION` for RNA).

### Incremental processing

The pipeline detects already-processed subjects by checking for existing output files. Subjects with all expected outputs are skipped automatically, allowing you to add new subjects to the samplesheet and re-run without reprocessing the entire cohort.

---

## Prerequisites

### Software

| Software | Version | Notes |
|---|---|---|
| [Nextflow](https://www.nextflow.io/) | >= 24.10.2 | `module load nextflow/25.10.2` on most HPC systems |
| [Apptainer](https://apptainer.org/) | >= 1.0 | Required for all container execution. Not Docker. |
| [htslib](http://www.htslib.org/) | Any recent | Provides `bcftools` and `tabix` |

### Container images

Five Apptainer containers are used. They are pulled automatically on first run and cached in `containers/`:

| Container | Purpose |
|---|---|
| `sdp-r:4.5.1` | R scripts (CNV, SV, expression, aggregation) |
| `sdp-python:3.12` | Python scripts |
| `vcf2maf + ensembl-vep` | VCF-to-MAF conversion |
| `pcgr:2.2.5` | Germline variant annotation (PCGR/CPSR) |
| `sdp-sigprofiler:1.1.3` | Mutational signature fitting (SigProfilerAssignment) |

> **GHCR authentication**: Containers hosted at `ghcr.io/crchum-citadel` require membership in the `crchum-citadel` GitHub organization and a personal access token:
> ```bash
> apptainer registry login --username <github_username> oras://ghcr.io
> ```

### Reference data

These must be pre-staged before running, especially on HPC compute nodes that lack internet access.

| Reference | Parameter | Description |
|---|---|---|
| GRCh38 FASTA | `genome_reference` | Human reference genome with `.fai` index. Required for vcf2maf. |
| VEP cache | `vep_data` | Ensembl VEP cache directory (version 113). If empty, the pipeline attempts to download a test cache. |
| PCGR reference | `pcgr_data` | PCGR reference data directory. If empty, the pipeline attempts to download from PCGR servers. |
| Ensembl annotations (CNV/SV) | `ensembl_annotations` | BioMart TSV (Ensembl 113) with `entrez_ncbi_id` column. Used for CNV gene mapping and SV breakpoint annotation. |
| Ensembl annotations (expression) | `ensembl_annotations_expr` | BioMart TSV (Ensembl 110) for Isofox Ensembl-to-Entrez ID mapping. |
| COSMIC fusion data | `cosmic_data` | COSMIC Fusion TSV (e.g., `Cosmic_Fusion_v103_GRCh38.tsv`) for ML SV processing. |
| ChimerKB data | `chimer_data` | ChimerDB KnownFusions Excel file. Bundled as `assets/test_data/annotations/ChimerKB4.xlsx`. |
| COSMIC v3.6 signatures | `cosmic_reference` | Bundled at `assets/COSMIC_Human_v3.6.zip`. No action needed unless overriding. |
| Signature metadata | `sbs_metadata`, `dbs_metadata`, `id_metadata` | Bundled in `assets/`. No action needed unless overriding. |

> **Offline HPC nodes**: Set `NXF_OFFLINE=true` and pre-pull all containers on a login node before submitting jobs:
> ```bash
> apptainer pull --dir containers/ oras://ghcr.io/crchum-citadel/sdp-r:4.5.1
> ```

---

## Preparing input data

### Genomic samplesheet

The genomic samplesheet is a CSV with one row per subject.

| Column | Required | Description |
|---|---|---|
| `group` | Yes | Study/cohort identifier (no spaces) |
| `subject_id` | Yes | Patient/subject identifier |
| `sample_id` | Yes | Sample base name used in output filenames |
| `folder` | Yes | Absolute path to the oncoanalyser output directory for this subject |

**Example:**

```csv
group,subject_id,sample_id,folder
COHORT1,PATIENT1,SAMPLE1,/data/oncoanalyser_output/COHORT1/PATIENT1
COHORT1,PATIENT2,SAMPLE2,/data/oncoanalyser_output/COHORT1/PATIENT2
COHORT2,PATIENT3,SAMPLE3,/data/oncoanalyser_output/COHORT2/PATIENT3
```

#### Expected folder layout

Each subject's `folder` must contain oncoanalyser output organized as follows. **All modalities are optional** — missing or empty files are silently skipped.

```
<folder>/
├── pave/
│   ├── <subject_id>-T.pave.somatic.vcf.gz          # Somatic DNA mutations
│   └── <subject_id>-T.pave.germline.vcf.gz         # Germline DNA mutations
├── sage_append/
│   └── somatic/
│       └── <subject_id>-T.sage.append.vcf.gz       # RNA-supported somatic mutations
├── esvee/
│   └── <subject_id>-T.esvee.somatic.vcf.gz         # Somatic structural variants
├── purple/
│   ├── <subject_id>-T.purple.cnv.somatic.tsv       # Somatic copy-number segments
│   └── <subject_id>-T.purple.cnv.gene.tsv          # Gene-level copy numbers
├── isofox/
│   ├── <subject_id>-T-RNA.isf.gene_data.csv        # Gene expression (Isofox)
│   └── <subject_id>-T-RNA.isf.pass_fusions.csv     # RNA gene fusions (Isofox)
└── sigs/
    └── <subject_id>-T.sig.snv_counts.csv           # Trinucleotide SNV counts (96-channel)
```

> **File naming**: The pipeline resolves files using `subject_id` (from the samplesheet), not `sample_id`. For example, if `subject_id=PATIENT1`, the pipeline looks for `PATIENT1-T.pave.somatic.vcf.gz`.

### Clinical samplesheet

The clinical samplesheet is a CSV pointing to individual clinical data files.

| Column | Required | Description |
|---|---|---|
| `group_id` | No | Group/study identifier |
| `filetype` | Yes | Data type (see table below) |
| `filepath` | Yes | Path to clinical CSV file |
| `extraction_date` | Yes | Date the data was extracted (YYYY-MM-DD) |
| `info` | No | Additional information |

**Supported filetypes:**

| Filetype | Description |
|---|---|
| `patient` | Patient demographics |
| `diagnosis` | Diagnosis data |
| `treatment` | Treatment records |
| `surgeries` | Surgical procedures |
| `systemic_treatment` | Systemic treatments (chemotherapy, immunotherapy, etc.) |
| `specimen` | Specimen/biobank data |
| `radio_therapy` | Radiotherapy records |

**Example:**

```csv
group_id,filetype,filepath,extraction_date,info
COHORT1,patient,/data/clinical/01_patient.csv,2025-05-20,
COHORT1,diagnosis,/data/clinical/02_diagnosis.csv,2025-05-20,
COHORT1,treatment,/data/clinical/03_treatment.csv,2025-05-20,
COHORT1,surgeries,/data/clinical/04_surgeries.csv,2025-05-20,
COHORT1,systemic_treatment,/data/clinical/05_systemic_treatment.csv,2025-05-20,
COHORT1,specimen,/data/clinical/06_specimen.csv,2025-05-20,
COHORT1,radio_therapy,/data/clinical/08_radio_therapy.csv,2025-05-20,
```

### ID linking file (optional)

A TSV file mapping genomic sample IDs to clinical patient IDs. Used when genomic and clinical data use different identifier schemes.

```
genomic_id	clinical_id
SAMPLE1	PATIENT_001
SAMPLE2	PATIENT_002
```

Pass via `--id_linking_file path/to/linking.tsv`.

---

## Configuration reference

All parameters are set in `nextflow.config` or overridden on the command line with `--param_name value`.

### Core parameters

| Parameter | Required | Default | Description |
|---|---|---|---|
| `mode` | Yes | `genomic` | Pipeline mode: `genomic` or `clinical` |
| `project_name` | Yes | — | Study name written to `meta_study.txt` |
| `project_description` | Yes | — | Study description written to `meta_study.txt` |
| `outdir` | No | `output_<project>_<date>` | Output directory |

### Genomic input parameters

| Parameter | Required | Default | Description |
|---|---|---|---|
| `genomic_samplesheet` | Yes | — | Path to genomic samplesheet CSV |
| `ensembl_annotations` | Yes | — | BioMart TSV (Ensembl 113) for CNV/SV gene annotation |
| `ensembl_annotations_expr` | Yes | — | BioMart TSV (Ensembl 110) for expression Entrez mapping |
| `genome_reference` | Yes | — | GRCh38 reference FASTA (for vcf2maf) |
| `vep_data` | No | (empty) | VEP cache directory. Empty = auto-download test cache |
| `vep_params` | No | `--cache-version 113 --ncbi-build=GRCh38` | Additional VEP CLI arguments |
| `pcgr_data` | No | (empty) | PCGR reference data directory. Empty = auto-download |
| `cosmic_data` | No | — | COSMIC Fusion TSV for ML SV processing |
| `chimer_data` | No | bundled | ChimerKB4.xlsx for ML fusion annotation |
| `cosmic_reference` | No | bundled | COSMIC v3.6 signature reference ZIP |
| `sbs_metadata` | No | bundled | SBS signature metadata TSV |
| `dbs_metadata` | No | bundled | DBS signature metadata TSV |
| `id_metadata` | No | bundled | ID (indel) signature metadata TSV |

### Clinical input parameters

| Parameter | Required | Default | Description |
|---|---|---|---|
| `clinical_samplesheet` | Yes | — | Path to clinical samplesheet CSV |
| `id_linking_file` | No | (empty) | TSV mapping genomic to clinical IDs |

### Container parameters

| Parameter | Default | Description |
|---|---|---|
| `container_r` | ORAS or local SIF | R 4.5.1 container |
| `container_python` | ORAS or local SIF | Python 3.12 container |
| `container_vcf2maf` | Wave container | vcf2maf + Ensembl-VEP |
| `container_pcgr` | ORAS or local SIF | PCGR 2.2.5 |
| `container_sigprofiler` | ORAS or local SIF | SigProfiler 1.1.3 |

### Notification parameters

| Parameter | Default | Description |
|---|---|---|
| `email` | — | Email address for pipeline completion notifications |
| `email_on_fail` | — | Email address for failure notifications |
| `plaintext_email` | `false` | Send plain-text emails instead of HTML |
| `hook_url` | `null` | Webhook URL for notifications (e.g., Slack) |

---

## Running the pipeline

### Profiles

Always include the `apptainer` profile. Combine with other profiles as needed.

| Profile | Purpose |
|---|---|
| `apptainer` | **Required.** Enables Apptainer containers, sets cache directory to `containers/`. |
| `slurm` | Submits jobs to SLURM. Default account: `def-chasse`. Edit `clusterOptions` in `nextflow.config` for your allocation. |
| `debug` | Dumps process hashes, prints hostname, disables work directory cleanup. |
| `gpu` | Adds `--nv` flag to Apptainer for GPU passthrough. |
| `local` | Uses local SIF file paths instead of ORAS URIs (for environments without GHCR access). |
| `test` | Uses bundled test data with reduced resource limits. Good for verifying installation. |

### Genomic mode

```bash
module load nextflow/25.10.2 apptainer htslib

nextflow run main.nf \
    -profile slurm,apptainer \
    --mode genomic \
    --project_name "MY_STUDY" \
    --project_description "Description of my study" \
    --genomic_samplesheet samplesheet.csv \
    --ensembl_annotations /path/to/biomart_grch38_ensembl_113.tsv \
    --ensembl_annotations_expr /path/to/biomart_grch38_ensembl_110.tsv \
    --genome_reference /path/to/Homo_sapiens_assembly38.fasta \
    --vep_data /path/to/vep_cache \
    --pcgr_data /path/to/pcgr_ref_data
```

### Clinical mode

```bash
nextflow run main.nf \
    -profile slurm,apptainer \
    --mode clinical \
    --project_name "MY_STUDY" \
    --project_description "Description of my study" \
    --clinical_samplesheet clinical_samplesheet.csv
```

### Using a params file

Instead of passing every parameter on the command line, create a YAML file:

```yaml
# params.yaml
mode: "genomic"
project_name: "MY_STUDY"
project_description: "Description of my study"
genomic_samplesheet: "/data/samplesheet.csv"
ensembl_annotations: "/path/to/biomart_grch38_ensembl_113.tsv"
ensembl_annotations_expr: "/path/to/biomart_grch38_ensembl_110.tsv"
genome_reference: "/path/to/Homo_sapiens_assembly38.fasta"
vep_data: "/path/to/vep_cache"
pcgr_data: "/path/to/pcgr_ref_data"
cosmic_data: "/path/to/Cosmic_Fusion_v103_GRCh38.tsv"
outdir: "output_my_study"
email: "user@example.com"
```

Then run:

```bash
nextflow run main.nf -profile slurm,apptainer -params-file params.yaml
```

> **Warning**: Use `-params-file` (not `-c`) for parameters. Custom config files via `-c` should only be used for process resource tuning and infrastructure settings.

### Local execution (no SLURM)

```bash
nextflow run main.nf -profile apptainer \
    --mode genomic \
    --genomic_samplesheet samplesheet.csv
```

### Running via Nextflow container

If Nextflow is not installed on your system, run it inside an Apptainer container:

```bash
apptainer pull --dir containers/ oras://ghcr.io/crchum-citadel/sdp-nextflow:25.04.7

apptainer exec containers/sdp-nextflow_25.04.7.sif \
    nextflow run main.nf -profile apptainer \
    --mode genomic \
    --genomic_samplesheet samplesheet.csv
```

### Resuming a failed run

Nextflow caches intermediate results. If a run fails, fix the issue and resume:

```bash
nextflow run main.nf -profile slurm,apptainer -params-file params.yaml -resume
```

Nextflow reuses cached results from completed steps and only re-runs what failed. Additionally, the pipeline's incremental processing skips subjects that already have all expected output files in the output directory.

### Running tests

Verify your installation with the bundled test data:

```bash
# Run all locally-testable tests
./nf-test test tests/clinical.nf.test tests/subworkflows/ tests/modules/ --profile test,apptainer

# Run a single test
./nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer
```

See the [README](../README.md#testing) for more test commands.

---

## Troubleshooting

### Container pull failures

**Symptom**: `FATAL: Unable to pull from ORAS` or authentication errors.

**Fix**: Log in to GHCR on the node that pulls containers:

```bash
apptainer registry login --username <github_username> oras://ghcr.io
```

For offline compute nodes, pre-pull containers on a login node:

```bash
apptainer pull --dir containers/ oras://ghcr.io/crchum-citadel/sdp-r:4.5.1
# Repeat for all 5 containers
```

### VEP/PCGR download failures

**Symptom**: VEP or PCGR steps fail trying to download reference data.

**Fix**: Pre-stage VEP cache and PCGR data on a login node with internet access, then set `vep_data` and `pcgr_data` to the staged paths. Compute nodes on most HPC clusters have no internet.

### SLURM account errors

**Symptom**: `sbatch: error: invalid account` or similar.

**Fix**: Edit `nextflow.config` and change the `clusterOptions` in the `slurm` profile to your allocation:

```groovy
clusterOptions = '--account=your-allocation'
```

### "All subjects in samplesheet already processed"

**Symptom**: Pipeline exits with error saying all subjects are already processed.

**Cause**: The incremental processing logic found all expected output files for every subject.

**Fix**: Either (a) delete the output directory to reprocess, (b) remove already-processed subjects from the samplesheet, or (c) add new subjects that haven't been processed yet.

### Expression merge requires 2+ samples

**Symptom**: Warning `Found 1 TPM file(s) for group X. Need at least 2 files to merge. Skipping merge step.`

**Cause**: cBioPortal expression files require at least 2 samples per group to produce a valid matrix.

**Fix**: Add more subjects to the group, or accept that the expression file will not be generated for single-sample groups.

### Memory or time errors

The pipeline uses process labels to set resource limits (see `conf/base.config`):

| Label | CPUs | Memory | Time |
|---|---|---|---|
| `process_single` | 1 | 6 GB | 4 h |
| `process_low` | 2 | 12 GB | 4 h |
| `process_medium` | 6 | 36 GB | 8 h |
| `process_high` | 12 | 72 GB | 16 h |
| `process_high_memory` | — | 200 GB | — |

Processes with the `error_retry` label automatically retry up to 3 times. To increase resources for a specific process, add a custom config:

```groovy
// custom.config
process {
    withName: 'VCF2MAF' {
        memory = '100.GB'
        time   = '24.h'
    }
}
```

Then pass it with `-c custom.config`.

### Missing input files

**Symptom**: Warnings like `File not found for PATIENT1 (mutation (PAVE somatic)): /path/to/file`.

**Cause**: The expected file does not exist in the subject's oncoanalyser output folder.

**Fix**: Check the folder layout matches the [expected structure](#expected-folder-layout). Missing files are skipped — the pipeline processes whatever modalities are available.

### Nextflow memory

If Nextflow itself runs out of memory, set this in your shell profile:

```bash
export NXF_OPTS='-Xms1g -Xmx4g'
```

### Running in the background

Use `screen`, `tmux`, or submit the Nextflow process itself as a SLURM job:

```bash
# Using screen
screen -S nf-run
nextflow run main.nf -profile slurm,apptainer -params-file params.yaml
# Ctrl+A, D to detach; screen -r nf-run to reattach
```
