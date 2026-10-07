# crchum-citadel/ampliseq-cbioportal

## Introduction

**crchum-citadel/ampliseq-cbioportal** formats ampliseq genomic data (VCFs + structural variant TSVs) into cBioPortal-compatible files: mutations (MAF), copy number alterations, structural variants, segmentation data, and clinical data.

Requires Nextflow >= 26.04.4 and Apptainer.

## Usage

### 1. Prepare input files

#### Generating a samplesheet

A helper script can automatically generate the samplesheet from a directory of sample folders:

```bash
python3 bin/generate_samplesheet.py /path/to/sample_data -o samplesheet.csv
```

The input directory can contain either cohort subdirectories (each with per-sample folders) or sample folders directly. Each sample folder must contain `*-basespace-pisces.final.vcf.gz` and `analysis_*_export.tsv`. The cohort subdirectory name becomes the `group` column, the `sample_id` is extracted from the VCF filename prefix (before `-basespace-`), and the `subject_id` is derived from the sample*id (part before the first `*`). Samples missing required files are skipped with warnings.

#### Samplesheet format

**Samplesheet** (`samplesheet.csv`):

```csv
group,subject_id,sample_id,folder_location
cohort_A,PATIENT_001,SAMPLE_001,path/to/samples/SAMPLE_001
```

Each sample folder must contain:

- `analysis_*_export.tsv` — structural variant / CNA export
- `*-basespace-pisces.final.vcf.gz` — compressed VCF for mutation calling
- `*-basespace-cnv.final.vcf` — CNV VCF for segmentation (`CN` FORMAT field required)

**Linking file** (`linking_file.txt`, tab-separated) — maps anonymized → real IDs:

```
sample_id	deanon_sample_id	deanon_patient_id
SAMPLE_001	PATIENT_001	PATIENT_001
```

One patient may have multiple rows (one per sample). `deanon_patient_id` is used to filter `data_clinical_patient.txt` to only patients whose samples are in the samplesheet.

**Patient file** (tab-separated): `patient_id`, `moh_id`, `age`, `sex`, `os_status`, `os_months`, `smoking_history`. `moh_id` is carried through to the `MOH_ID` column of `data_clinical_patient.txt`.

**Sample file** (tab-separated): `sample_id`, `patient_id`, `cancer_type`, `cancer_type_detailed`, `sample_type`, `primary_tumor_site`, `metastatic_tumor_site`, `tumor_purity`

All listed columns are required. Other columns are ignored. `Unknown` and empty values become `NA`. See `assets/patient_file.txt` and `assets/sample_file.txt` for examples.

> The `sample_id` column in the sample file must use the **deanonymized** (real) sample IDs — the same values that appear in the `deanon_sample_id` column of the linking file. Outputs are scoped to the samplesheet: only samples present in the samplesheet appear in `data_clinical_sample.txt`, `data_clinical_patient.txt`, and `case_lists/`, even if the sample file and patient file contain additional entries.

### 2. Run the pipeline

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results/ \
  --patient_file patient_file.txt \
  --sample_file sample_file.txt \
  --linking_file linking_file.txt \
  --mafsmith_data /path/to/mafsmith_home/ \
  --study_id my_study
```

VCF → MAF uses [mafsmith](https://github.com/nf-osi/mafsmith) + fastVEP (GRCh37, PASS records only), run from the `mafsmith-fastvep` container (see [Containers](#containers)).

`--mafsmith_data` is optional. When it is omitted, empty or `null`, the GRCh37 reference bundle (Ensembl 113) is downloaded once into `assets/mafsmith/` and reused by later runs. A pre-staged bundle must contain `GRCh37/reference.fa` and `GRCh37/genes.gff3.gz`. `MAFSMITH` fails if either file is missing or empty.

> A `-stub-run` with `--skip_vcf2maf false` stores an **empty** `assets/mafsmith/`. Delete it afterwards, or later real runs will reuse it and fail on the missing reference files.

Skip VCF → MAF conversion if MAFs already exist:

```bash
nextflow run main.nf ... --skip_vcf2maf true
```

Filter mode applies QC thresholds using the sample's `analysis_*_export.tsv`. By default it is off and everything is kept:
```bash
nextflow run main.nf ... --filter_tsv_variants true
```

| Data | Kept in filter mode | Params (default) |
|---|---|---|
| Mutations | Position overlaps a TSV row with `Depth` ≥ 250 and `VAF` > 0.03 (`chr` prefix ignored) | `--mutation_min_depth` (250), `--mutation_min_vaf` (0.03) |
| CNAs | `Confidence` = HIGH and raw `Copy Number` ≥ 6, so all deletions are dropped | `--cna_confidence` (HIGH), `--cna_min_copy_number` (6) |
| Fusions | Taken from the TSV, even when a `*-star-fusion.final.vcf` exists. `Supporting Reads` summed per gene pair ≥ 1000 | `--sv_min_supporting_reads` (1000) |

`data_seg.txt` is not filtered. Filter mode runs per sample, so samples already in `--outdir` are not re-filtered. Use a fresh `--outdir` when you turn it on or change a threshold.

Boolean flags accept `true`/`false` on the command line. Nextflow 26 passes CLI values as strings, so the pipeline converts them explicitly.

Resume a previous run:

```bash
nextflow run main.nf ... -resume
```

### Containers

| Label | Param | Default | Definition |
|---|---|---|---|
| `python` | `--python_sif` | local `.sif` | `containers/python-ampliseq.def` |
| `mafsmith` | `--mafsmith_container` | `oras://ghcr.io/crchum-citadel/mafsmith-fastvep:0.1.0-0.4.0` | `containers/mafsmith-fastvep_v0.1.0-0.4.0.def` |

The ghcr images are private; log in first with `apptainer registry login`. To build locally instead, see [`containers/README.md`](containers/README.md).

### Testing

```bash
./nf-test test tests/modules --profile=+apptainer   # all module tests (MAFSMITH pulls the mafsmith image)
nextflow run main.nf -profile test,apptainer        # end-to-end on assets/ test data (skip_vcf2maf=true)
```

### Incremental runs

When new samples are added to the samplesheet, re-running with the same `--outdir` automatically skips samples already processed. Only new samples go through per-sample steps; the merge and downstream steps re-run over all samples.

```bash
# Initial run
nextflow run main.nf --input samplesheet_v1.csv --outdir results/ ...

# Later — add new samples to the samplesheet and re-run
nextflow run main.nf --input samplesheet_v2.csv --outdir results/ ...
# → existing samples skipped, new samples processed, outputs updated
```

A sample is skipped when all four files exist under `results/samples/{sample_id}/`:
`{sample_id}_sv.txt`, `{sample_id}_cna.txt`, `{sample_id}_seg.txt`, `{sample_id}_mutations.txt`

### Running standalone scripts (without Nextflow)

The full transformation can also be run via:

```bash
# Edit hardcoded paths at the top of the script first
bash bin/run_pipeline.sh
```

Or run individual scripts from within the output directory:

```bash
cd /path/to/output

python3 /path/to/bin/format_tsv_to_sv.py       <analysis_export.tsv> <SAMPLE_ID>
python3 /path/to/bin/format_cna.py       <analysis_export.tsv> <SAMPLE_ID>
python3 /path/to/bin/format_mutations.py data_mutations.txt    <linking_file>
python3 /path/to/bin/format_sv.py        data_sv.txt           <linking_file>
python3 /path/to/bin/format_cna_deanon.py data_cna.txt         <linking_file>
python3 /path/to/bin/vcf_to_seg.py        <cnv.vcf>             <SAMPLE_ID>
python3 /path/to/bin/seg_deanon.py        data_seg.txt          <linking_file>
python3 /path/to/bin/clinical_patients_format.py <patient_file> <linking_file>
python3 /path/to/bin/clinical_sample_format.py   <sample_file> <linking_file>
```
