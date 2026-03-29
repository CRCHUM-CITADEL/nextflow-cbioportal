# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow pipeline (`crchum-citadel/ampliseq-cbioportal`) that formats ampliseq genomic data into cBioPortal-compatible format. Built from nf-core/tools template v3.5.1. The workflow DAG is implemented in `workflows/ampliseq-cbioportal.nf`, split into three subworkflows under `subworkflows/local/`, with each process in its own file under `modules/local/`. The core transformation logic lives in standalone scripts in `bin/`, which are also directly runnable outside Nextflow via `bin/run_pipeline.sh`.

## Running the Pipeline

```bash
# Run with required parameters
nextflow run main.nf --input samplesheet.csv --outdir results/ \
  --patient_file patient_file.txt \
  --sample_file sample_file.txt \
  --linking_file linking_file.txt \
  --vcf2maf_container community.wave.seqera.io/library/vcf2maf_ensembl-vep:... \
  --vep_data /path/to/vep_data/ \
  --study_id optilab_study

# Skip VCF → MAF if MAFs already exist
nextflow run main.nf ... --skip_vcf2maf true

# Pass all mutations through without TSV-coordinate filtering
nextflow run main.nf ... --filter_tsv_variants false

# Resume a previous run
nextflow run main.nf ... -resume
```

Requires Nextflow >= 25.04.0. All container paths are Apptainer `.sif` images.

## Generating a Samplesheet

Use `bin/generate_samplesheet.py` to auto-build `samplesheet.csv` from a data directory:

```bash
python3 bin/generate_samplesheet.py <input_dir> [--output samplesheet.csv]
```

`input_dir` should contain cohort subdirectories, each holding per-sample folders. The script derives `group` from the cohort directory name, `sample_id` from the `*-basespace-pisces.final.vcf.gz` filename prefix, and `subject_id` from the part of `sample_id` before the first `_`. If `input_dir` contains sample folders directly (no cohort subdirectories), it is treated as a single cohort named after the directory. Folders missing either required file are skipped with a warning.

## Running the Standalone Transformation Scripts

As an alternative to Nextflow, the full transformation can be run via `bin/run_pipeline.sh`. The script contains **hardcoded cluster paths** (`BASE_DIR`, `SCRIPTS_DIR`, `DATA_DIR`, `APPTAINER_SIF`, `LINKING_FILE`) that must be updated for each environment before running.

```bash
# Edit paths at the top of the script, then:
bash bin/run_pipeline.sh
```

All Python scripts write output relative to `os.getcwd()`, so they must be run from the target output directory. The orchestrator does this via `(cd "$OUT_DIR" && python3.12 script.py ...)`.

```bash
# Run individual scripts from within the output directory
cd /path/to/output

python3 /path/to/bin/format_tsv.py    <analysis_export.tsv> <SAMPLE_ID>   # appends to data_sv.txt
python3 /path/to/bin/format_cna.py    <analysis_export.tsv> <SAMPLE_ID>   # appends to data_cna.txt
python3 /path/to/bin/format_mutations.py data_mutations.txt <linking_file> # deanonymizes in-place
python3 /path/to/bin/format_sv.py       data_sv.txt         <linking_file> # deanonymizes in-place
python3 /path/to/bin/format_cna_deanon.py data_cna.txt      <linking_file> # deanonymizes in-place
python3 /path/to/bin/clinical_patients_format.py <patient_file>            # writes data_clinical_patient.txt
python3 /path/to/bin/clinical_sample_format.py   <sample_file>             # writes data_clinical_sample.txt
python3 /path/to/bin/format_meta.py <study_id> [out_dir]                 # writes all meta_*.txt files
```

## Input File Formats

**Pipeline samplesheet** (`--input`, `test_data/samplesheet.csv`):
```
group,subject_id,sample_id,folder_location
cohort_A,PATIENT_001,SAMPLE_001,test_data/samples/SAMPLE_001
```
Note: `assets/samplesheet.csv` uses the nf-core FASTQ schema (`sample,fastq_1,fastq_2`) — this is a template artifact, not the actual input format.

**Per-sample folder** must contain:
- `analysis_*_export.tsv` — tab-separated with columns: `Chr`, `Start`, `End`, `Variant Type`, `Variant Subtype`, `Genes`, `Breakend Genes`, `Supporting Reads`, `Copy Number`
- `*-basespace-pisces.final.vcf.gz` — compressed VCF; filename prefix becomes the `SAMPLE_ID`

**Linking file** (`linking_file.txt`) — tab-separated, maps anonymized → real IDs:
```
sample_id	deanon_sample_id
SAMPLE_001	PATIENT_001
```

**Patient file** — tab-separated: `patient_id`, `age`, `sex`, `os_status` (0/1), `os_months`, `smoking_history`

**Sample file** — tab-separated: `num_id`, `sample_id`, `patient_id`, `cancer_type`, `cancer_type_detailed`, `sample_type`, `tumor_site`, `tumor_purity`

## Architecture

```
main.nf                                  # Entry point; reads --input CSV, calls AMPLISEQ_CBIOPORTAL workflow
workflows/
  ampliseq-cbioportal.nf                 # Full workflow DAG: per-sample → merge → deanon → clinical → meta
subworkflows/local/
  per_sample_format/main.nf             # Per-sample: SV, CNA, mutations
  merge_deanon/main.nf                  # Merge per-sample files and deanonymize
  study_metadata/main.nf               # Clinical files, case lists, meta files
modules/local/
  format_sv/main.nf                     # Extracts FUSION rows → _sv.txt
  format_cna/main.nf                    # Extracts DUPLICATION/DELETION rows → _cna.txt
  vcf_to_maf/main.nf                   # Runs vcf2maf via Apptainer container
  stub_maf/main.nf                      # Emits empty MAF header (skip_vcf2maf=true)
  filter_mutations/main.nf             # Filters MAF rows by TSV coordinates (filter_tsv_variants=true)
  passthrough_mutations/main.nf        # Copies MAF through without filtering (filter_tsv_variants=false)
  merge_sv/main.nf                      # Concatenates per-sample _sv.txt files
  merge_cna/main.nf                     # Concatenates per-sample _cna.txt files
  merge_mutations/main.nf              # Concatenates per-sample _mutations.txt files
  deanon_mutations/main.nf             # Deanonymizes Tumor_Sample_Barcode in data_mutations.txt
  deanon_sv/main.nf                    # Deanonymizes Sample_Id in data_sv.txt
  deanon_cna/main.nf                   # Deanonymizes Sample_Id in data_cna.txt
  clinical_patients/main.nf            # Formats patient file → data_clinical_patient.txt
  clinical_samples/main.nf             # Formats sample file → data_clinical_sample.txt
  write_case_lists/main.nf             # Writes case_lists/ from linking file
  write_meta/main.nf                   # Writes all cBioPortal meta_*.txt files
bin/
  generate_samplesheet.py              # Auto-builds samplesheet.csv from a data directory
  run_pipeline.sh                       # Orchestrates full transformation outside Nextflow (hardcoded cluster paths)
  format_tsv.py                         # Extracts FUSION rows → data_sv.txt
  format_cna.py                         # Extracts DUPLICATION/DELETION rows → data_cna.txt (long format)
  format_mutations.py                  # Deanonymizes Tumor_Sample_Barcode in data_mutations.txt
  format_sv.py                          # Deanonymizes Sample_Id in data_sv.txt
  format_cna_deanon.py                 # Deanonymizes Sample_Id in data_cna.txt
  clinical_patients_format.py          # Formats patient file → data_clinical_patient.txt
  clinical_sample_format.py            # Formats sample file → data_clinical_sample.txt
  format_meta.py                        # Writes all cBioPortal meta_*.txt files for a study
assets/
  samplesheet.csv                       # Test samplesheet pointing to assets/samples/
  linking_file.txt / patient_file.txt / sample_file.txt  # Test clinical files
  samples/SAMPLE_00{1,2}/              # Each has analysis_*_export.tsv and *.vcf.gz
  schema_input.json                    # JSON schema for samplesheet validation (nf-core template)
nextflow.config                        # Process defaults, profiles, all pipeline params
nextflow_schema.json                   # Parameter schema for --help and validation
```

## Data Flow

**Per-sample:**
1. `analysis_*_export.tsv` → FORMAT_SV → `_sv.txt` (FUSION rows)
2. `analysis_*_export.tsv` → FORMAT_CNA → `_cna.txt` (DUPLICATION/DELETION rows)
3. VCF → VCF_TO_MAF (vcf2maf, VEP v113, GRCh37/hg19) → MAF → FILTER_MUTATIONS (filter by TSV coordinates, default) or PASSTHROUGH_MUTATIONS (skip filtering) → `_mutations.txt`

**Downstream:**
4. MERGE_SV / MERGE_CNA / MERGE_MUTATIONS collect per-sample files into merged files
5. DEANON_MUTATIONS / DEANON_SV / DEANON_CNA replace anonymized IDs using linking file → `data_mutations.txt`, `data_sv.txt`, `data_cna.txt`
6. CLINICAL_PATIENTS + CLINICAL_SAMPLES write cBioPortal 5-line-header format files
7. WRITE_CASE_LISTS generates `case_lists/` from deanonymized IDs in linking file
8. WRITE_META writes cBioPortal study meta files (`meta_study.txt`, `meta_mutations.txt`, `meta_sv.txt`, `meta_cna.txt`, `meta_clinical_patient.txt`, `meta_clinical_sample.txt`)

Output files: `data_mutations.txt`, `data_sv.txt`, `data_cna.txt`, `data_clinical_patient.txt`, `data_clinical_sample.txt`, `case_lists/`

## Key Implementation Notes

- `format_tsv.py` and `format_cna.py` both read the same `analysis_*_export.tsv` but filter on different `Variant Subtype` values: `FUSION` (SVs) vs `DUPLICATION`/`DELETION` (CNAs)
- CNA copy number → cBioPortal value mapping (TSV mode): CN=2 is normal and dropped; non-2 values map to `0→-2, 1→-1, 3→1, ≥4→2`
- `data_cna.txt` is written in long format (Hugo_Symbol, Sample_Id, Value); `meta_cna.txt` declares `datatype: DISCRETE_LONG` so cBioPortal accepts this format directly — no pivot needed
- The vcf2maf Apptainer container mounts `vep_data` as `/home/jbellavance/` inside the container
- `clinical_sample_format.py` reads only the first 8 columns of the sample file and drops `num_id` and `tumor_purity`
- All deanon scripts warn to stderr on unmatched IDs and leave them unchanged
- `filter_tsv_variants` controls mutation filtering only: when `true` (default) mutations are filtered to TSV coordinates; when `false` all MAF mutations pass through unfiltered

## Development Status

First stable release. The full Nextflow pipeline is functional end-to-end. Each process is in its own file under `modules/local/`, organized into three subworkflows. GitHub CI/CD, nf-test, MultiQC, and nf-core documentation were intentionally skipped from the nf-core template.
