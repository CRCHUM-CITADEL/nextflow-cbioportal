# CRCHUM-CITADEL/nextflow-cbioportal

## Clone

To run this pipeline, you first need to git clone this repository and enter the directory:

```
git clone https://github.com/CRCHUM-CITADEL/nextflow-cbioportal.git && cd
```

> [!NOTE]
> For all containers we are using Apptainer because of it's compatibility with HPC environments.
> If you're not working from an HPC environment, find out how to install it here: https://apptainer.org/docs/admin/main/installation.html

## Change nextflow.config file

You will need to change parameters in the nextflow config in order to point to certain files. These options are found in the `params` dict in nextflow.config. Parameters are mandatory unless specified otherwise.

| Field                    | Description                                                                                                                    |
| ------------------------ | ------------------------------------------------------------------------------------------------------------------------------ |
| mode                     | Pipeline run mode. Options: `clinical`, `genomic`                                                                              |
| genomic_samplesheet      | Input samplesheet for genomic pipeline. See section below.                                                                     |
| ensembl_annotations_expr | Ensembl annotation .tsv file for expression subworkflow (tested with ensembl 110 with biomart)                                 |
| ensembl_annotations      | Ensembl annotation .tsv file. (tested with 113 with biomart)                                                                   |
| vep_cache                | Cache folder of downloaded ensembl vep release.                                                                                |
| vep_params               | Parameters for VEP usage as described here:`<br>` https://github.com/Ensembl/ensembl-vep?tab=readme-ov-file#options (optional) |
| pcgr_data                | Folder of pcgr reference data (uncompressed)                                                                                   |
| genome_reference         | Location of GRCh38 reference fasta file.                                                                                       |
| container_pcgr           | Location of PCGR apptainer image (remote or local)                                                                             |
| container_python         | Location of Python apptainer image (remote or local)                                                                           |
| container_r              | Location of R apptainer image (remote or local)                                                                                |
| container_vcf2maf        | Location of nf-core vcf2maf module container`<br>` apptainer image (remote or local) (optional)                                |
| clinical_samplesheet     | Input samplesheet for clinical pipeline. See section below.                                                                    |
| id_linking_file          | ID linking file generated from genomic pipeline.                                                                               |

## Samplesheet

You will need to create a samplesheet for this pipeline, which can differ between modes.

### Generating a samplesheet

A helper script can automatically generate the genomic samplesheet from DRAGEN output:

```bash
python3 create_samplesheet.py --input_dir /path/to/dragen_output
```

This expects a directory containing `dna/` and `rna/` subfolders with standard DRAGEN tumor pair output. It produces `samplesheet.csv` with columns: `group_id`, `subject_id`, `sample_id`, `sample_type`, `sequence_data`, `info`, `filepath`. Three rows are created per subject (DT somatic DNA, DN germinal DNA, RT somatic RNA).

To run the script's built-in validation tests:

```bash
python3 create_samplesheet.py --test
```

### Mode = 'genomic'

The samplesheet format is heavily based on `<a href="https://github.com/nf-core/oncoanalyser">` oncoanalyser's nf-core pipeline `</a>`. See below for exact specfications:

#### Genomic Input Schema

The genomic input file must be a CSV file where each object contains the following fields:

| Column Name     | Type    | Required | Pattern                      | Options               | Description                                                                        |
| --------------- | ------- | -------- | ---------------------------- | --------------------- | ---------------------------------------------------------------------------------- |
| `group_id`      | string  | No       | `^\S+$` (no spaces)          | -                     | Group identifier                                                                   |
| `subject_id`    | string  | **Yes**  | `^(?:\d+\|\S+)$` (no spaces) | -                     | Subject identifier                                                                 |
| `sample_id`     | integer | **Yes**  | `^\d+$` (numeric only)       | -                     | Sample identifier                                                                  |
| `sample_type`   | string  | **Yes**  | -                            | `somatic`, `germinal` | Type of sample                                                                     |
| `sequence_data` | string  | **Yes**  | -                            | `dna`, `rna`          | Type of sequence data                                                              |
| `info`          | string  | No       | -                            | -                     | Additional information                                                             |
| `filepath`      | string  | **Yes**  | -                            | -                     | Path to DRAGEN output folder containing sample germinal or tumoral data (DN/DT/RT) |

> [!NOTE]
> Fields marked as **Required** must be present in each object
> All string fields cannot contain spaces unless otherwise noted

### mode = 'clinical'

#### Clinical Input Schema

The clinical input file must be a CSV where each object contains the following fields:

| Column Name       | Type   | Required | Pattern             | Options                                                                                                   | Description                                            |
| ----------------- | ------ | -------- | ------------------- | --------------------------------------------------------------------------------------------------------- | ------------------------------------------------------ |
| `group_id`        | string | No       | `^\S+$` (no spaces) | -                                                                                                         | Group identifier                                       |
| `filetype`        | string | **Yes**  | -                   | `patient`, `diagnosis`, `treatment`, `surgeries` `<br>` `systemic_treatment`, `specimen`, `radio_therapy` | File type category                                     |
| `filepath`        | string | **Yes**  | `^\S+\.csv)$`       | -                                                                                                         | Path to clincial data file (must end with `.csv`)      |
| `extraction_date` | string | **Yes**  | `^\S+$`             | -                                                                                                         | Date of which the data was extracted from the database |
| `info`            | string | No       | -                   | -                                                                                                         | Additional information                                 |

> [!NOTE]
> Fields marked as **Required** must be present in each object
> All string fields cannot contain spaces unless otherwise noted
> The `filepath` must point to a valid file with one of the accepted extensions

TBD

## Running the pipeline for non-HPC environments (i.e. without SLURM)

### Pull nextflow container image

In order to run nextflow with the exact software used to build the pipeline, pull the container image hosted on CITADEL's organisational Github (or if not in a member in crchum-citadel GitHub, ask for the location of the locally-stored .sif file.)

```
apptainer pull --dir containers/ oras://ghcr.io/crchum-citadel/sdp-nextflow:25.04.7
```

> [!NOTE]
> in order to pull from the github repository, you need to be have your credentials stored in your environment, and be a member of the crchum-citadel GitHub.
> You will need to first create a `<a href="https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens" target="_blank">`PAT token with appropriate permissions`</a>`
> To authenticate:
> `apptainer registry login --username <username> oras://ghcr.io`
> You will be prompted for your PAT token.

### Run nextflow via apptainer

Structure your command as such:

```
apptainer exec containers/sdp-nextflow_<version>.sif nextflow <command>
```

To get the version:

```
apptainer exec containers/sdp-nextflow_v25.10.2.sif nextflow -v
```

To run the pipeline test (using minimal data):

```
apptainer exec containers/sdp-nextflow_v25.10.2.sif nextflow run main.nf -profile test,apptainer
```

## Running the pipeline in an HPC environment (i.e. with SLURM)

### Module load nextflow

```
module load nextflow/25.10.2 apptainer/1.3.5 htslib/1.22.1
```

### Run nextflow

```
nextflow run main.nf -profile apptainer,slurm
```

> [!IMPORTANT]
> Always run the pipeline with atleast the apptainer profile option, or else it won't work. (E.g. -profile apptainer)

## Nf-core functionality

This pipeline was build from a nf-core template. You may want to use some nf-core CLI functionality.

According to the documentation, you can run directly from the container:

```
apptainer exec \
    --bind $(pwd):$(pwd) \
    --pwd $(pwd) \
    docker://nfcore/tools:3.3.2 \
    nf-core pipelines list
```

> [!NOTE]
> You will likely want to create an alias in your ~/.bashrc file. In this file, include:
> `alias nf-core="apptainer exec --bind $(pwd):$(pwd) --pwd $(pwd) docker://nfcore/tools:3.3.2 nf-core"`
> (Don't forget to source and change versions when there's an update!)

> [!NOTE]
> If you're looking to update the container registry, you need to:
>
> 1. Build the image
>    `apptainer build <new_sif_file>.sif <def_file>.def`
> 2. Push the image:
>    `apptainer push <new_sif_file>.sif oras://ghcr.io/crchum-citadel/<name>:<version>`

## For developers:

> [!IMPORTANT]
> Always make changes in the dev branch of the repository. To merge changes to the main branch, create a pull request and make sure all tests pass.
> Alternatively, create an issue in the GitHub repository.

### Running tests

The test suite includes 8 tests: 1 pipeline-level clinical test and 7 subworkflow tests.

```bash
# Run all tests (~5 min)
./nf-test test tests/subworkflows/ tests/clinical.nf.test --profile test,apptainer

# Run individual subworkflow test
./nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer

# Update snapshots after intentional output change
./nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer --update-snapshot
```

| Test | Description | Time |
|------|-------------|------|
| `clinical.nf.test` | Pipeline-level clinical mode | ~15s |
| `clinical_aggregate.nf.test` | Clinical file aggregation | ~10s |
| `genomic_cnv.nf.test` | CNV SEG + long files from DRAGEN VCF | ~10s |
| `genomic_sv.nf.test` | SV files from fusion_candidates.final | ~10s |
| `genomic_expression.nf.test` | TPM from Sailfish quant.genes.sf | ~10s |
| `genomic_aggregate_output.nf.test` | Merged genomic files + meta/case-lists | ~10s |
| `genomic_ml.nf.test` | ML tables (stub-run) | ~10s |
| `genomic_mutations.nf.test` | VEP + PCGR + vcf2maf chain | ~200s |

### Other functionalities:

To run pre-commit (to check linting before pull request):

```
apptainer exec containers/nextflow-citadel_v25.10.2.sif pre-commit run .
```
