# CRCHUM-CITADEL/nextflow-cbioportal (oncoanalyser): Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- Timeline generation from ARGO clinical CSVs (`gen_timeline.R`): produces a single combined `data_timeline.txt` with all event types (surgery, treatment, status, specimen, lab_test)
- MOHCCN mapping tables (`assets/mohccn_*` CSVs) for translating plain-text values to ICD-O / ontology codes in clinical and timeline output
- New params: `mohccn_primary_site_map`, `mohccn_specimen_tissue_source_map`, `mohccn_treatment_intent_map`
- MOHCCN code integration in `clin_format.R` (primary site, tissue source, treatment intent)
- Specimen SAMPLE_TYPE derivation from nearest follow-up disease status (Primary / Recurrence / Metastasis)
- Biomarker lab test pivoting to cBioPortal long format (PSA, CEA, CA125, ER/PR/HER2, HPV)
- `GENERATE_TIMELINE` module in `clinical_aggregate` subworkflow
- DNA vs RNA SV distinction (`ecb509b`): DNA_Support/RNA_Support flags properly set per source
- `combine_cbioportal_outputs.py`: support for `data_timeline.txt` (union-of-columns append), `cancer_type.txt`, per-subject folders, and glob-discovered `meta_*.txt` / `case_lists/*.txt` so new outputs are never silently dropped
- `combine_cbioportal_outputs.py`: `--strict` flag, and a warning listing any input file no merge strategy handles
- `tests/test_combine_cbioportal_outputs.py` pytest suite, plus a pytest job in the linting workflow

### Fixed

- `<study>.tar.gz` in `both` mode was genomic-only: `PACKAGE_CBIOPORTAL` now runs from `main.nf` after `CLINICAL`, so the archive also carries the clinical files, `data_timeline.txt`, `meta_timeline.txt` and `util_linking_file.txt`
- `data_mutational_signatures_counts_ID.txt` was published to the output directory but missing from `<study>.tar.gz`
- `combine_cbioportal_outputs.py`: clinical and timeline merges no longer misalign rows when the two batches carry different column sets
- `combine_cbioportal_outputs.py`: a genuine expression value of `0` is no longer replaced by the fill value in wide-matrix merges
- `combine_cbioportal_outputs.py`: row-append no longer glues two rows together when the first file lacks a trailing newline (`util_linking_file.txt`)
- `combine_cbioportal_outputs.py`: in-place merge keeps a backup until the swap succeeds instead of deleting the target first
- `tests/.nftignore` excludes `*.tar.gz`: tar embeds per-file mtimes, so the archive md5 differed on every run
- Timeline START_DATE defaults to 0 (diagnosis day) when source date is missing
- Clinical sample deduplication: better merging to remove duplicate sample-level rows
- Removed `analyte_type` from sample-level clinical output rows
- Status data frame handling in clinical output

## v1.0.0dev - [date]

Initial release, created with the [nf-core](https://nf-co.re/) template.

### Added

- Genomic mode: mutations (vcf2maf), CNV (PURPLE), SV (ESVEE + Isofox fusions), expression (Isofox), mutational signatures (SBS/DBS/ID)
- Clinical mode: ARGO/ICGC-ARGO CSV processing into cBioPortal patient/sample attribute files
- Both mode: automatic linking of genomic sample IDs to clinical patient IDs
- Incremental processing: skip already-processed subjects on re-runs
- ML feature table generation (CNV, expression, mutation, SV)
- Samplesheet generators for genomic and clinical input
- cBioPortal output packaging (.tar.gz)
