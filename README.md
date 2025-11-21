# CRCHUM-CITADEL/nextflow-cbioportal (English)

(Traduction en francais suit)

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

| Field                | Description                                                                                                                   |
| -------------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| mode                 | Pipeline run mode. Options : ['clinical', 'genomic']                                                                          |
| genomic_samplesheet  | Input samplesheet for genomic pipeline. See section below.                                                                    |
| ensembl_annotations_expr  | Ensembl annotation .tsv file for expression subworkflow (tested with ensembl 110 with biomart)                                                                                                 |
| ensembl_annotations  | Ensembl annotation .tsv file. (tested with 113 with biomart)                                                                                                 |
| vep_cache            | Cache folder of downloaded ensembl vep release.                                                                               |
| vep_params           | Parameters for VEP usage as described here: <br> https://github.com/Ensembl/ensembl-vep?tab=readme-ov-file#options (optional) |
| pcgr_data            | Folder of pcgr reference data (uncompressed)                                                                                  |
| genome_reference     | Location of GRCh38 reference fasta file.                                                                                      |
| container_pcgr       | Location of PCGR apptainer image (remote or local)                                                                            |
| container_python     | Location of Python apptainer image (remote or local)                                                                          |
| container_r          | Location of R apptainer image (remote or local)                                                                               |
| container_vcf2maf    | Location of nf-core vcf2maf module container <br> apptainer image (remote or local) (optional)                                |
| clinical_samplesheet | Input samplesheet for clinical pipeline. See section below.                                                                   |
| id_linking_file      | ID linking file generated from genomic pipeline.                                                                              |

## Samplesheet

You will need to create a samplesheet for this pipeline, which can differ between modes.

### Mode = 'genomic'

The samplesheet format is heavily based on <a href="https://github.com/nf-core/oncoanalyser"> oncoanalyser's nf-core pipeline </a>. See below for exact specfications:

#### Genomic Input Schema

The genomic input file must be a CSV file where each object contains the following fields:

| Column Name     | Type    | Required | Pattern                         | Options                                    | Description                                                             |
| --------------- | ------- | -------- | ------------------------------- | ------------------------------------------ | ----------------------------------------------------------------------- |
| `group_id`      | string  | No       | `^\S+$` (no spaces)             | -                                          | Group identifier                                                        |
| `subject_id`    | string  | **Yes**  | `^(?:\d+\|\S+)$` (no spaces)    | -                                          | Subject identifier                                                      |
| `sample_id`     | integer | **Yes**  | `^\d+$` (numeric only)          | -                                          | Sample identifier                                                       |
| `sample_type`   | string  | **Yes**  | -                               | `somatic`, `germinal`                      | Type of sample                                                          |
| `sequence_data` | string  | **Yes**  | -                               | `dna`, `rna`                               | Type of sequence data                                                   |
| `filetype`      | string  | **Yes**  | -                               | `cnv`, `sv`, `expression`, `hard_filtered` | File type category                                                      |
| `info`          | string  | No       | -                               | -                                          | Additional information                                                  |
| `filepath`      | string  | **Yes**  | `^\S+\.(?:vcf\.gz\|final\|sf)$` | -                                          | Path to genomic data file (must end with `.vcf.gz`, `.final`, or `.sf`) |

> [!NOTE]
> Fields marked as **Required** must be present in each object
> All string fields cannot contain spaces unless otherwise noted
> The `filepath` must point to a valid file with one of the accepted extensions

### mode = 'clinical'

#### Clinical Input Schema

The clinical input file must be a CSV where each object contains the following fields:

| Column Name       | Type   | Required | Pattern             | Options                                                                                                 | Description                                            |
| ----------------- | ------ | -------- | ------------------- | ------------------------------------------------------------------------------------------------------- | ------------------------------------------------------ |
| `group_id`        | string | No       | `^\S+$` (no spaces) | -                                                                                                       | Group identifier                                       |
| `filetype`        | string | **Yes**  | -                   | `patient`, `diagnosis`, `treatment`, `surgeries` <br> `systemic_treatment`, `specimen`, `radio_therapy` | File type category                                     |
| `filepath`        | string | **Yes**  | `^\S+\.csv)$`       | -                                                                                                       | Path to clincial data file (must end with `.csv`)      |
| `extraction_date` | string | **Yes**  | `^\S+$`             | -                                                                                                       | Date of which the data was extracted from the database |
| `info`            | string | No       | -                   | -                                                                                                       | Additional information                                 |

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
> You will need to first create a <a href="https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens" target="_blank">PAT token with appropriate permissions</a>  
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
apptainer exec containers/sdp-nextflow_v25.04.7.sif nextflow -v
```

To run the pipeline test (using minimal data):

```
apptainer exec containers/sdp-nextflow_v25.04.7.sif nextflow run main.nf -profile test,apptainer
```

## Running the pipeline in an HPC environment (i.e. with SLURM)

### Module load nextflow

```
module load nextflow/25.04.6 apptainer
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
apptainer exec containers/nextflow-citadel_v25.04.7.sif nf-test test --profile apptainer
```

To run pre-commit (to check linting before pull request):

```
apptainer exec containers/nextflow-citadel_v25.04.7.sif pre-commit run .
```

# CRCHUM-CITADEL/nextflow-cbioportal (Francais)

Documentation de démarrage (Traduction en français)

## Cloner

Pour exécuter ce pipeline, vous devez d'abord cloner ce dépôt git et entrer dans le répertoire :

```
git clone https://github.com/CRCHUM-CITADEL/nextflow-cbioportal.git && cd
```

> [!NOTE]
> Pour tous les conteneurs, nous utilisons Apptainer en raison de sa compatibilité avec les environnements HPC.
> Si vous ne travaillez pas dans un environment HPC, trouvez comment l'installer ici : https://apptainer.org/docs/admin/main/installation.html

## Modifier le fichier nextflow.config

Vous devrez modifier les paramètres dans le fichier de configuration nextflow afin de pointer vers certains fichiers. Ces options se trouvent dans le dictionnaire `params` dans nextflow.config. Obligatoire sauf indication contraire.

| Champ                | Description                                                                                                                                |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------ |
| mode                 | Mode d'exécution du pipeline. Options : ['clinical', 'genomic']                                                                            |
| genomic_samplesheet  | Feuille d'échantillons génomique en entrée. Voir la section ci-dessous.                                                                    |
| gencode_annotations  | Fichier d'annotations Gencode .gtf.                                                                                                        |
| ensembl_annotations  | Fichier d'annotations Ensembl .tsv.                                                                                                        |
| vep_cache            | Dossier du cache de la version ensembl vep téléchargée.                                                                                    |
| vep_params           | Paramètres pour l'utilisation de VEP comme décrit ici : <br> https://github.com/Ensembl/ensembl-vep?tab=readme-ov-file#options (optionnel) |
| pcgr_data            | Dossier des données de référence pcgr (décompressées)                                                                                      |
| genome_reference     | Emplacement du fichier fasta de référence GRCh38.                                                                                          |
| container_pcgr       | Emplacement de l'image apptainer PCGR (distant ou local)                                                                                   |
| container_python     | Emplacement de l'image apptainer Python (distant ou local)                                                                                 |
| container_r          | Emplacement de l'image apptainer R (distant ou local)                                                                                      |
| container_vcf2maf    | Emplacement de l'image apptainer du module <br> nf-core vcf2maf (distant ou local) (optionnel)                                             |
| clinical_samplesheet | Feuille d'échantillons clinique en entrée. Voir la section ci-dessous.                                                                     |
| id_linking_file      | Fichier de linkage d'échantillion, créé par pipeline génomique                                                                             |

## Feuille d'échantillons

Vous devrez créer une feuille d'échantillons pour ce pipeline, qui peut différer selon les modes.

### Mode = 'genomic'

Le format de la feuille d'échantillons est fortement basé sur <a href="https://github.com/nf-core/oncoanalyser"> le pipeline nf-core d'oncoanalyser </a>. Voir ci-dessous pour les spécifications exactes :

#### Schéma d'entrée génomique

Le fichier d'entrée génomique doit être un tableau JSON où chaque objet contient les champs suivants :

| Nom de colonne  | Type   | Requis  | Motif                            | Options                                    | Description                                                                                      |
| --------------- | ------ | ------- | -------------------------------- | ------------------------------------------ | ------------------------------------------------------------------------------------------------ |
| `group_id`      | chaîne | Non     | `^\S+$` (pas d'espaces)          | -                                          | Identifiant de groupe                                                                            |
| `subject_id`    | chaîne | **Oui** | `^(?:\d+\|\S+)$` (pas d'espaces) | -                                          | Identifiant du sujet                                                                             |
| `sample_id`     | entier | **Oui** | `^\d+$` (numérique uniquement)   | -                                          | Identifiant d'échantillon                                                                        |
| `sample_type`   | chaîne | **Oui** | -                                | `somatic`, `germinal`                      | Type d'échantillon                                                                               |
| `sequence_data` | chaîne | **Oui** | -                                | `dna`, `rna`                               | Type de données de séquençage                                                                    |
| `filetype`      | chaîne | **Oui** | -                                | `cnv`, `sv`, `expression`, `hard_filtered` | Catégorie de type de fichier                                                                     |
| `info`          | chaîne | Non     | -                                | -                                          | Informations supplémentaires                                                                     |
| `filepath`      | chaîne | **Oui** | `^\S+\.(?:vcf\.gz\|final\|sf)$`  | -                                          | Chemin vers le fichier de données génomiques (doit se terminer par `.vcf.gz`, `.final` ou `.sf`) |

> [!NOTE]
> Les champs marqués comme **Requis** doivent être présents dans chaque objet
> Tous les champs de type chaîne ne peuvent pas contenir d'espaces sauf indication contraire
> Le `filepath` doit pointer vers un fichier valide avec l'une des extensions acceptées

### mode = 'clinical'

#### Schéma d’entrée clinique

Le fichier d’entrée clinique doit être un fichier CSV où chaque objet contient les champs suivants :

| Nom de colonne    | Type   | Requis  | Motif                  | Options                                                                                                 | Description                                                               |
| ----------------- | ------ | ------- | ---------------------- | ------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- |
| `group_id`        | chaîne | Non     | `^\S+$` (sans espaces) | -                                                                                                       | Identifiant de groupe                                                     |
| `filetype`        | chaîne | **Oui** | -                      | `patient`, `diagnosis`, `treatment`, `surgeries` <br> `systemic_treatment`, `specimen`, `radio_therapy` | Catégorie du type de fichier                                              |
| `filepath`        | chaîne | **Oui** | `^\S+\.csv)$`          | -                                                                                                       | Chemin vers le fichier de données cliniques (doit se terminer par `.csv`) |
| `extraction_date` | chaîne | **Oui** | `^\S+$`                | -                                                                                                       | Date à laquelle les données ont été extraites de la base de données       |
| `info`            | chaîne | Non     | -                      | -                                                                                                       | Informations supplémentaires                                              |

> [!NOTE]
> Les champs marqués comme **Requis** doivent être présents dans chaque objet.  
> Les champs de type chaîne ne peuvent pas contenir d’espaces sauf indication contraire.  
> Le `filepath` doit pointer vers un fichier valide avec une des extensions acceptées.

## Pour les environnements non-HPC (c'est-à-dire sans SLURM)

### Télécharger l'image conteneur nextflow

Afin d'exécuter nextflow avec le logiciel exact utilisé pour construire le pipeline, téléchargez l'image conteneur hébergée sur le GitHub organisationnel de CITADEL (ou si vous n'êtes pas membre de crchum-citadel GitHub, demandez l'emplacement du fichier .sif stocké localement.)

```
apptainer pull --dir containers/ oras://ghcr.io/crchum-citadel/sdp-nextflow:25.04.7
```

> [!NOTE]
> Pour télécharger depuis le dépôt GitHub, vous devez avoir vos identifiants stockés dans votre environnement et être membre du GitHub crchum-citadel.
> Vous devrez d'abord créer un <a href="https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens" target="_blank">jeton PAT avec les permissions appropriées</a>  
> Pour vous authentifier :
> `apptainer registry login --username <nom_utilisateur> oras://ghcr.io`
> Vous serez invité à entrer votre jeton PAT.

### Exécuter nextflow via apptainer

Structurez votre commande comme suit :

```
apptainer exec containers/sdp-nextflow_<version>.sif nextflow <commande>
```

Pour obtenir la version :

```
apptainer exec containers/sdp-nextflow_v25.04.7.sif nextflow -v
```

Pour exécuter le test du pipeline (en utilisant des données minimales) :

```
apptainer exec containers/sdp-nextflow_v25.04.7.sif nextflow run main.nf -profile test,apptainer
```

## Dans un environnement HPC (c'est-à-dire avec SLURM)

### Charger le module nextflow

```
module load nextflow/25.04.6 apptainer
```

### Exécuter nextflow

```
nextflow run main.nf -profile apptainer,slurm
```

> [!IMPORTANT]
> Exécutez toujours le pipeline avec au moins l'option de profil apptainer, sinon il ne fonctionnera pas. (Par ex. -profile apptainer)

## Fonctionnalité nf-core

Ce pipeline a été construit à partir d'un modèle nf-core. Vous voudrez peut-être utiliser certaines fonctionnalités CLI de nf-core.

Selon la documentation, vous pouvez exécuter directement depuis le conteneur :

```
apptainer exec \
    --bind $(pwd):$(pwd) \
    --pwd $(pwd) \
    docker://nfcore/tools:3.3.2 \
    nf-core pipelines list
```

> [!NOTE]
> Vous voudrez probablement créer un alias dans votre fichier ~/.bashrc. Dans ce fichier, incluez :
> `alias nf-core="apptainer exec --bind $(pwd):$(pwd) --pwd $(pwd) docker://nfcore/tools:3.3.2 nf-core"`
> (N'oubliez pas de sourcer et de changer les versions lors d'une mise à jour !)

> [!NOTE]
> Si vous cherchez à mettre à jour le registre de conteneurs, vous devez :
>
> 1. Construire l'image
>    `apptainer build <nouveau_fichier_sif>.sif <fichier_def>.def`
> 2. Pousser l'image :
>    `apptainer push <nouveau_fichier_sif>.sif oras://ghcr.io/crchum-citadel/<nom>:<version>`

## Pour les développeurs :

> [!IMPORTANT]
> Effectuez toujours les modifications dans la branche dev du dépôt. Pour fusionner les modifications vers la branche main, créez une demande de tirage (pull request) et assurez-vous que tous les tests réussissent.
> Alternativement, créez un problème (issue) dans le dépôt GitHub.

### Autres fonctionnalités :

Pour exécuter nf-test :

```
apptainer exec containers/nextflow-citadel_v25.04.7.sif nf-test test --profile apptainer
```

Pour exécuter pre-commit (pour vérifier le linting avant la demande de tirage) :

```
apptainer exec containers/nextflow-citadel_v25.04.7.sif pre-commit run .
```
