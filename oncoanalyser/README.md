# oncoanalyser → cBioPortal

Nextflow DSL2 pipeline that converts [oncoanalyser](https://github.com/nf-core/oncoanalyser) WGS/WTS output and clinical CSVs into [cBioPortal](https://www.cbioportal.org/)-ready files. Runs on HPC (SLURM + Apptainer — no Docker, no sudo).

> [!NOTE]
> All containers use Apptainer for compatibility with HPC environments.
> To install Apptainer: https://apptainer.org/docs/admin/main/installation.html

---

## Clone

```bash
git clone https://github.com/CRCHUM-CITADEL/nextflow-cbioportal.git
cd nextflow-cbioportal/oncoanalyser
```

---

## Configuration

Edit `nextflow.config` to point to your data. Parameters are mandatory unless marked optional.

| Parameter | Description |
|---|---|
| `mode` | Run mode: `'genomic'` or `'clinical'` |
| `genomic_samplesheet` | CSV samplesheet for genomic mode (one row per subject) |
| `clinical_samplesheet` | CSV samplesheet for clinical mode |
| `id_linking_file` | TSV linking genomic sample IDs to clinical IDs |
| `ensembl_annotations` | BioMart TSV for gene annotation (GRCh38, must include `entrez_ncbi_id`) |
| `vep_data` | Path to pre-staged VEP cache directory (empty → downloads test cache) |
| `genome_reference` | GRCh38 reference FASTA |
| `container_r` | R Apptainer image URI |
| `container_python` | Python Apptainer image URI |
| `container_vcf2maf` | vcf2maf Apptainer image URI (optional — mutations mode only) |
| `container_pcgr` | PCGR Apptainer image URI (optional — mutations mode only) |
| `cosmic_data` | COSMIC/ChimerDB fusion data file for ML step |
| `outdir` | Output directory (default: `output/`) |

---

## Samplesheets

### Genomic mode (`--mode genomic`)

One row per subject. The `folder` column points to that subject's oncoanalyser output directory.

| Column | Required | Description |
|---|---|---|
| `group` | Yes | Study/cohort identifier (no spaces) |
| `subject_id` | Yes | Patient/subject identifier |
| `sample_id` | Yes | Sample base name used in output filenames |
| `folder` | Yes | Path to the oncoanalyser output directory for this subject |

```csv
group,subject_id,sample_id,folder
COHORT1,PATIENT1,SAMPLE1,/data/oncoanalyser_output/COHORT1/PATIENT1
COHORT1,PATIENT2,SAMPLE2,/data/oncoanalyser_output/COHORT1/PATIENT2
```

#### Expected folder layout

The pipeline resolves files from each subject's `folder` using this layout:

```
<folder>/
├── pave/
│   ├── <subject>-T.pave.somatic.vcf.gz           ← somatic DNA mutations
│   └── <subject>-T.pave.germline.vcf.gz          ← germline DNA mutations
├── sage_append/
│   └── somatic/
│       └── <subject>-T.sage.append.vcf.gz        ← RNA-supported somatic mutations
├── esvee/
│   └── <subject>-T.esvee.unfiltered.vcf.gz       ← somatic structural variants (tumor only)
├── purple/
│   ├── <subject>-T.purple.cnv.somatic.tsv        ← somatic copy-number segments
│   └── <subject>-T.purple.cnv.gene.tsv           ← gene-level copy numbers
└── isofox/
    ├── <subject>-T-RNA.isf.gene_data.csv         ← gene expression (Isofox)
    └── <subject>-T-RNA.isf.pass_fusions.csv      ← RNA-seq gene fusions (Isofox)
```

> [!NOTE]
> Missing or empty files are silently skipped — each modality is optional. DNA SVs and RNA fusions are both written to the group-level `data_sv.txt`; they share the same column schema and are distinguished by the `DNA_Support` / `RNA_Support` flags and `Class` (`DELETION`/`DUPLICATION`/`INVERSION`/`TRANSLOCATION` for DNA, `FUSION` for RNA).

### Clinical mode (`--mode clinical`)

| Column | Required | Description |
|---|---|---|
| `group_id` | No | Group identifier |
| `filetype` | Yes | One of: `patient`, `diagnosis`, `treatment`, `surgeries`, `systemic_treatment`, `specimen`, `radio_therapy` |
| `filepath` | Yes | Path to clinical CSV file |
| `extraction_date` | Yes | Date the data was extracted from source |
| `info` | No | Additional information |

---

## Running the pipeline

### HPC (SLURM + Apptainer)

```bash
module load nextflow/25.10.2 apptainer htslib

nextflow run main.nf -profile slurm,apptainer \
    --mode genomic \
    --genomic_samplesheet samplesheet.csv
```

> [!IMPORTANT]
> Always include the `apptainer` profile: `-profile apptainer` (or `slurm,apptainer`)

### Local (Apptainer, no SLURM)

```bash
nextflow run main.nf -profile apptainer \
    --mode genomic \
    --genomic_samplesheet samplesheet.csv
```

### Via Nextflow container

```bash
apptainer pull --dir containers/ oras://ghcr.io/crchum-citadel/sdp-nextflow:25.04.7

apptainer exec containers/sdp-nextflow_25.04.7.sif \
    nextflow run main.nf -profile apptainer --mode genomic --genomic_samplesheet samplesheet.csv
```

> [!NOTE]
> Pulling from GHCR requires membership in the `crchum-citadel` GitHub organisation and a PAT token:
> ```bash
> apptainer registry login --username <username> oras://ghcr.io
> ```

---

## Containers

Images are pulled automatically on first run and cached in `containers/`. To pull manually or build from source:

```bash
# Pull from GHCR (requires auth)
apptainer pull --dir containers/ oras://ghcr.io/crchum-citadel/sdp-r:4.5.1

# Build locally from .def file (no sudo needed with --fakeroot)
apptainer build --fakeroot containers/ghcr.io-crchum-citadel-sdp-r-4.5.1.sif containers/r_v4.5.1.def

# Create .img symlink so Nextflow's ORAS cache lookup succeeds
ln -sf ghcr.io-crchum-citadel-sdp-r-4.5.1.sif containers/ghcr.io-crchum-citadel-sdp-r-4.5.1.img
```

To push a new image to the registry:

```bash
apptainer push <image>.sif oras://ghcr.io/crchum-citadel/<name>:<version>
```

---

## Testing

Tests use [nf-test](https://www.nf-test.com/) (binary at `./nf-test`):

```bash
# Run all locally-testable tests
./nf-test test tests/clinical.nf.test tests/subworkflows/ tests/modules/ --profile test,apptainer

# Run a single subworkflow test
./nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer

# Run module tests
./nf-test test tests/modules/isofox_fusion_to_cbioportal.nf.test --profile test,apptainer

# Update snapshots after intentional output change
./nf-test test tests/subworkflows/genomic_cnv.nf.test --profile test,apptainer --update-snapshot
```

Test fixtures live in `assets/test_data/`. To regenerate them from scratch:

```bash
bash bin/create_test_data.sh
```

CI runs all tests on GitHub Actions with sharding (up to 7 parallel jobs). The `genomic_mutations` test requires PCGR and vcf2maf containers and is CI-only.

---

## Development

> [!IMPORTANT]
> Always work in the `dev` branch. Merge to `main` via pull request after tests pass.

### Linting

```bash
pre-commit run --all-files
```

### nf-core tools

```bash
# Run via container (no local install needed)
apptainer exec \
    --bind $(pwd):$(pwd) --pwd $(pwd) \
    docker://nfcore/tools:3.3.2 \
    nf-core pipelines list

# Or add an alias in ~/.bashrc:
# alias nf-core="apptainer exec --bind $(pwd):$(pwd) --pwd $(pwd) docker://nfcore/tools:3.3.2 nf-core"

nf-core modules install <module_name>
nf-core subworkflows install <subworkflow_name>
```
