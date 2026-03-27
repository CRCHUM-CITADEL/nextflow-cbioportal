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
| mode                     | Pipeline run mode. Options : ['clinical', 'genomic']                                                                           |
| genomic_samplesheet      | Input samplesheet for genomic pipeline (one row per subject). See section below.                                               |
| ensembl_annotations_expr | Ensembl annotation .tsv file for expression subworkflow (tested with ensembl 110 with biomart)                                 |
| ensembl_annotations      | Ensembl annotation .tsv file. (tested with 113 with biomart)                                                                   |
| vep_data                 | Cache folder of downloaded ensembl vep release.                                                                                |
| vep_params               | Parameters for VEP usage as described here:`<br>` https://github.com/Ensembl/ensembl-vep?tab=readme-ov-file#options (optional) |
| genome_reference         | Location of GRCh38 reference fasta file.                                                                                       |
| container_python         | Location of Python apptainer image (remote or local)                                                                           |
| container_r              | Location of R apptainer image (remote or local)                                                                                |
| container_vcf2maf        | Location of nf-core vcf2maf module container`<br>` apptainer image (remote or local) (optional)                                |
| clinical_samplesheet     | Input samplesheet for clinical pipeline. See section below.                                                                    |
| id_linking_file          | ID linking file generated from genomic pipeline.                                                                               |
| cosmic_data              | Cosmic fusion data .tsv                                                                                                         |

## Samplesheet

You will need to create a samplesheet for this pipeline, which can differ between modes.

### Mode = 'genomic'

The genomic samplesheet has **one row per subject**. The `folder` column points to the oncoanalyser output directory for that subject; all modality files are resolved relative to it.

#### Genomic Input Schema

| Column Name  | Type   | Required | Pattern                      | Description                                               |
| ------------ | ------ | -------- | ---------------------------- | --------------------------------------------------------- |
| `group`      | string | **Yes**  | `^\S+$` (no spaces)          | Study/cohort group identifier                             |
| `subject_id` | string | **Yes**  | `^(?:\d+\|\S+)$` (no spaces) | Patient/subject identifier                                |
| `sample_id`  | string | **Yes**  | `^\S+$` (no spaces)          | Sample base name used in output file names                |
| `folder`     | string | **Yes**  | -                            | Path to the oncoanalyser output directory for this subject |

> [!NOTE]
> All string fields cannot contain spaces.

Example:

```csv
group,subject_id,sample_id,folder
COHORT1,PATIENT1,SAMPLE1,/data/oncoanalyser_output/COHORT1/PATIENT1
COHORT1,PATIENT2,SAMPLE2,/data/oncoanalyser_output/COHORT1/PATIENT2
```

#### Expected folder layout

The pipeline resolves output files from each subject's `folder` using this layout:

```
<folder>/
├── sage/
│   ├── somatic/
│   │   └── <sample_id>-T.sage.somatic.vcf.gz     ← somatic DNA mutations
│   ├── germline/
│   │   └── <sample_id>-N.sage.germline.vcf.gz    ← germline DNA mutations
│   └── append/
│       └── <sample_id>-T.sage.append.vcf.gz      ← somatic RNA mutations
├── esvee/
│   └── caller/
│       └── <sample_id>-T.esvee.unfiltered.vcf.gz ← structural variants
├── purple/
│   ├── <sample_id>-T.purple.cnv.somatic.tsv
│   └── <sample_id>-T.purple.cnv.gene.tsv
├── <sample_id>-T.isf.fusions.csv                 ← RNA fusions (Isofox)
└── <sample_id>-T.isf.gene_data.csv               ← gene expression (Isofox)
```

> [!NOTE]
> Files that are absent or empty are skipped with a warning — each modality is optional.

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
module load nextflow/25.10.2 apptainer htslib
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

## For developpers:

> [!IMPORTANT]
> Always make changes in the dev branch of the repository. To merge changes to the main branch, create a pull request and make sure all tests pass.
> Alternatively, create an issue in the GitHub repository.

### Other functionalities:

To run nf-test:

```
apptainer exec containers/nextflow-citadel_v25.10.2.sif nf-test test --profile apptainer
```

To run pre-commit (to check linting before pull request):

```
apptainer exec containers/nextflow-citadel_v25.10.2.sif pre-commit run .
```
