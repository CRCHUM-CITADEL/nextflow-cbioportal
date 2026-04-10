# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This is a Nextflow pipeline (`crchum-citadel/ampliseq-cbioportal`) that formats ampliseq genomic data into cBioPortal-compatible format, including mutations (MAF), discrete copy number alterations, structural variants, segmentation data (SEG), and clinical files. Built from nf-core/tools template v3.5.1. The workflow DAG is implemented in `workflows/ampliseq-cbioportal.nf`, split into three subworkflows under `subworkflows/local/`, with each process in its own file under `modules/local/`. The core transformation logic lives in standalone scripts in `bin/`, which are also directly runnable outside Nextflow via `bin/run_pipeline.sh`.

## Running the Pipeline

```bash
# Run with required parameters
nextflow run main.nf --input samplesheet.csv --outdir results/ \
  --patient_file patient_file.txt \
  --sample_file sample_file.txt \
  --linking_file linking_file.txt \
  --vcf2maf_sif /path/to/vcf2maf_ensembl-vep*.sif \
  --vep_data /path/to/vep_data/ \
  --study_id optilab_study

# Skip VCF → MAF if MAFs already exist
nextflow run main.nf ... --skip_vcf2maf true

# Resume a previous run
nextflow run main.nf ... -resume
```

Requires Nextflow >= 25.04.0. All container paths are Apptainer `.sif` images.

## Incremental Runs

When a cohort grows (new samples added to the samplesheet), re-running with the same `--outdir` automatically skips samples whose per-sample outputs already exist. Only new samples are processed; the merge/deanon/clinical steps re-run over all samples combined.

```bash
# Initial run — processes SAMPLE_001, SAMPLE_002
nextflow run main.nf --input samplesheet_v1.csv --outdir results/ ...

# Incremental run — skips SAMPLE_001 and SAMPLE_002, processes SAMPLE_003 only
nextflow run main.nf --input samplesheet_v2.csv --outdir results/ ...
```

A sample is considered already processed when all four files exist under `{outdir}/samples/{sample_id}/`:
`{sample_id}_sv.txt`, `{sample_id}_cna.txt`, `{sample_id}_seg.txt`, `{sample_id}_mutations.txt`

Skipped samples are logged: `INFO: Skipping already-processed sample: SAMPLE_001`

**Important:** use the same `--outdir` across runs for incremental skipping to work. Per-sample outputs are published with `mode: 'copy'` (not symlinks), so they survive `work/` cleanup.

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
python3 /path/to/bin/vcf_to_seg.py     <cnv.vcf>            <SAMPLE_ID>   # appends to data_seg.txt
python3 /path/to/bin/seg_deanon.py     data_seg.txt         <linking_file> # deanonymizes in-place
python3 /path/to/bin/clinical_patients_format.py <patient_file> <linking_file>                 # writes data_clinical_patient.txt (filtered to samplesheet patients)
python3 /path/to/bin/clinical_sample_format.py   <sample_file> <linking_file>                 # writes data_clinical_sample.txt (filtered to samplesheet samples)
python3 /path/to/bin/format_meta.py <study_id> [out_dir]                 # writes all meta_*.txt files
```

## Input File Formats

**Pipeline samplesheet** (`--input`, `assets/samplesheet.csv`):
```
group,subject_id,sample_id,folder_location
cohort_A,PATIENT_001,SAMPLE_001,assets/samples/SAMPLE_001
```
`assets/samplesheet_sample1.csv` contains only SAMPLE_001 and is used as the first-run input in the incremental test.

**Per-sample folder** must contain:
- `analysis_*_export.tsv` — tab-separated with columns: `Chr`, `Start`, `End`, `Variant Type`, `Variant Subtype`, `Genes`, `Breakend Genes`, `Supporting Reads`, `Copy Number`
- `*-basespace-pisces.final.vcf.gz` — compressed VCF for mutation calling; filename prefix becomes the `SAMPLE_ID`
- `*-basespace-cnv.final.vcf` — uncompressed CNV VCF for segmentation; must have `CN` in FORMAT and `END` in INFO for non-point segments

**Linking file** (`linking_file.txt`) — tab-separated, maps anonymized → real IDs:
```
sample_id	deanon_sample_id	deanon_patient_id
SAMPLE_001	PATIENT_001	PATIENT_001
```
`deanon_sample_id` maps to the real sample ID used in clinical and data files. `deanon_patient_id` maps to the patient ID used in the patient clinical file. One patient may have multiple rows (multiple samples).

**Patient file** — tab-separated: `patient_id`, `age`, `sex`, `os_status` (0/1), `os_months`, `smoking_history`

**Sample file** — tab-separated: `num_id`, `sample_id`, `patient_id`, `cancer_type`, `cancer_type_detailed`, `sample_type`, `tumor_site`, `tumor_purity`. The `sample_id` column must use the **deanonymized** sample IDs (matching `deanon_sample_id` in the linking file), not the anonymized IDs used in the samplesheet.

## Architecture

```
main.nf                          # Entry point; reads --input CSV, calls AMPLISEQ_CBIOPORTAL workflow
workflows/
  ampliseq-cbioportal.nf         # Full workflow DAG: per-sample → merge → deanon → clinical → case lists → meta
modules/local/
  format_sv/main.nf                     # Extracts FUSION rows → _sv.txt; publishes to samples/{id}/
  format_cna/main.nf                    # Extracts DUPLICATION/DELETION rows → _cna.txt; publishes to samples/{id}/
  vcf_to_maf/main.nf                   # Runs vcf2maf via Apptainer container
  stub_maf/main.nf                      # Emits empty MAF header (skip_vcf2maf=true)
  filter_mutations/main.nf             # Filters MAF rows by TSV coordinates (filter_tsv_variants=true); publishes to samples/{id}/
  passthrough_mutations/main.nf        # Copies MAF through without filtering (filter_tsv_variants=false); publishes to samples/{id}/
  vcf_to_seg/main.nf                   # Converts *-basespace-cnv.final.vcf → per-sample _seg.txt; publishes to samples/{id}/
  filter_linking/main.nf               # Filters linking file to only samplesheet samples (by anonymized sample_id)
  merge_sv/main.nf                      # Concatenates per-sample _sv.txt files
  merge_cna/main.nf                     # Concatenates per-sample _cna.txt files
  merge_mutations/main.nf              # Concatenates per-sample _mutations.txt files
  merge_seg/main.nf                     # Concatenates per-sample _seg.txt files
  deanon_mutations/main.nf             # Deanonymizes Tumor_Sample_Barcode in data_mutations.txt
  deanon_sv/main.nf                    # Deanonymizes Sample_Id in data_sv.txt
  deanon_cna/main.nf                   # Deanonymizes Sample_Id in data_cna.txt
  deanon_seg/main.nf                   # Deanonymizes ID in data_seg.txt
  clinical_patients/main.nf            # Formats patient file → data_clinical_patient.txt (filtered to samplesheet patients)
  clinical_samples/main.nf             # Formats sample file → data_clinical_sample.txt (filtered to samplesheet samples)
  write_case_lists/main.nf             # Writes case_lists/ from filtered linking file (samplesheet samples only)
  write_meta/main.nf                   # Writes all cBioPortal meta_*.txt files
bin/
  generate_samplesheet.py              # Auto-builds samplesheet.csv from a data directory
  run_pipeline.sh                       # Orchestrates full transformation outside Nextflow (hardcoded cluster paths)
  test_incremental.sh                   # Two-phase bash test for incremental run behavior
  format_tsv.py                         # Extracts FUSION rows → data_sv.txt
  format_cna.py                         # Extracts DUPLICATION/DELETION rows → data_cna.txt (long format)
  format_mutations.py                  # Deanonymizes Tumor_Sample_Barcode in data_mutations.txt
  format_sv.py                          # Deanonymizes Sample_Id in data_sv.txt
  format_cna_deanon.py                 # Deanonymizes Sample_Id in data_cna.txt
  vcf_to_seg.py                         # Converts CNV VCF → data_seg.txt (PASS records; seg.mean = log2(CN/2))
  seg_deanon.py                         # Deanonymizes ID column in data_seg.txt
  clinical_patients_format.py          # Formats patient file → data_clinical_patient.txt
  clinical_sample_format.py            # Formats sample file → data_clinical_sample.txt
  format_meta.py                        # Writes all cBioPortal meta_*.txt files for a study
assets/
  schema_input.json              # JSON schema for samplesheet validation (nf-core template)
  samplesheet.csv                # Test samplesheet: SAMPLE_001 + SAMPLE_002
  samplesheet_sample1.csv        # Test samplesheet: SAMPLE_001 only (used in incremental test)
nextflow.config                  # Process defaults, profiles, all pipeline params
nextflow_schema.json             # Parameter schema for --help and validation
```

## Data Flow

**Per-sample:**
1. `analysis_*_export.tsv` → FORMAT_SV → `_sv.txt` (FUSION rows)
2. `analysis_*_export.tsv` → FORMAT_CNA → `_cna.txt` (DUPLICATION/DELETION rows)
3. VCF → VCF_TO_MAF (vcf2maf, VEP v113, GRCh37/hg19) → MAF → FILTER_MUTATIONS (filter by TSV coordinates, default) or PASSTHROUGH_MUTATIONS (skip filtering) → `_mutations.txt`
4. `*-basespace-cnv.final.vcf` → VCF_TO_SEG → `_seg.txt` (PASS records only; seg.mean = log2(CN/2))

**Downstream:**
5. Per-sample outputs are published to `{outdir}/samples/{sample_id}/` (`_sv.txt`, `_cna.txt`, `_seg.txt`, `_mutations.txt`); on re-runs samples with all four files present are skipped automatically
6. FILTER_LINKING filters the linking file to only rows whose anonymized `sample_id` appears in the samplesheet → `linking_filtered.txt`; this filtered file is used for all downstream steps
7. MERGE_SV / MERGE_CNA / MERGE_MUTATIONS / MERGE_SEG collect per-sample files (new + existing) into merged files
8. DEANON_MUTATIONS / DEANON_SV / DEANON_CNA / DEANON_SEG replace anonymized IDs using filtered linking file → `data_mutations.txt`, `data_sv.txt`, `data_cna.txt`, `data_seg.txt`
9. CLINICAL_SAMPLES writes `data_clinical_sample.txt` filtered to samples in the samplesheet (matched via `deanon_sample_id` in filtered linking)
10. CLINICAL_PATIENTS writes `data_clinical_patient.txt` filtered to patients whose `patient_id` appears in the filtered linking file's `deanon_patient_id` column
11. WRITE_CASE_LISTS generates `case_lists/` from the filtered linking file (samplesheet samples only)
12. WRITE_META writes cBioPortal study meta files (`meta_study.txt`, `meta_mutations.txt`, `meta_sv.txt`, `meta_cna.txt`, `meta_seg.txt`, `meta_clinical_patient.txt`, `meta_clinical_sample.txt`)

Output files: `data_mutations.txt`, `data_sv.txt`, `data_cna.txt`, `data_seg.txt`, `data_clinical_patient.txt`, `data_clinical_sample.txt`, `case_lists/`, `samples/` (per-sample cache)

## Key Implementation Notes

- `format_tsv.py` and `format_cna.py` both read the same `analysis_*_export.tsv` but filter on different `Variant Subtype` values: `FUSION` (SVs) vs `DUPLICATION`/`DELETION` (CNAs)
- CNA copy number → cBioPortal value mapping: `0→-2, 1→-1, 3→1, ≥4→2` (CN=2 is normal, yields `None` and is dropped)
- `data_cna.txt` is written in long format (Hugo_Symbol, Sample_Id, Value); `meta_cna.txt` declares `datatype: DISCRETE_LONG` so cBioPortal accepts this format directly — no pivot needed
- The vcf2maf Apptainer container mounts `vep_data` as `/home/jbellavance/` inside the container
- `clinical_sample_format.py <sample_file> <linking_file>` reads the first 8 columns of the sample file, filters rows to those whose `sample_id` is in the filtered linking file's `deanon_sample_id` column, and drops `num_id`; the `sample_id` column in the sample file must use the deanonymized (real) sample IDs — i.e. the same IDs that appear in `deanon_sample_id` of the linking file
- `clinical_patients_format.py <patient_file> <linking_file>` filters the patient file to only patients whose `patient_id` appears in the filtered linking file's `deanon_patient_id` column
- All deanon scripts warn to stderr on unmatched IDs and leave them unchanged
- `filter_tsv_variants` controls mutation filtering only: when `true` (default) mutations are filtered to TSV coordinates; when `false` all MAF mutations pass through unfiltered
- `vcf_to_seg.py` reads `*-basespace-cnv.final.vcf`, keeps only PASS records, parses `END` from INFO (defaults to POS for point variants), reads integer `CN` from the FORMAT/sample columns, and computes `seg.mean = log2(CN/2)`; CN=0 yields −3.0 as a homozygous-deletion sentinel; `num.mark` is always 1 since ampliseq VCFs carry no probe-count information
- `meta_seg.txt` declares `datatype: SEG` and `show_profile_in_analysis_tab: false`; cBioPortal uses SEG files for the copy-number segment viewer

## Development Status

The full Nextflow pipeline is functional end-to-end with incremental run support. All `bin/` scripts are implemented and used as Nextflow processes via `modules/local/`. Incremental behavior is tested via `bin/test_incremental.sh`. GitHub CI/CD, nf-test, MultiQC, and nf-core documentation were intentionally skipped from the nf-core template.
