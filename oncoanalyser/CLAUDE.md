# CLAUDE.md — oncoanalyser

Oncoanalyser WGS/WTS + clinical CSVs → cBioPortal. HPC-only (SLURM + Apptainer).

## Key Rules

- R scripts → `container_r`; Python → `container_python`; SigProfiler → `container_sigprofiler`
- `data_sv.txt` rows require Hugo symbols at both sites — filter unannotated rows
- SV classification (`gen_esvee_sv_to_cbioportal.R`): BND ALT strand → `(+,-)` DEL, `(-,+)` DUP, `(+,+)/(−,−)` INV, diff chr TRANSLOC. DNA SVs: `DNA_Support=Yes, RNA_Support=No`
- RNA fusions (`gen_isofox_fusion_to_cbioportal.R`): `Class=FUSION, DNA_Support=No, RNA_Support=Yes`. Both merge into `data_sv.txt`
- Isofox 3.0 `pass_fusions.tsv` packs both gene symbols into one `Name` column
  (`<up>_<down>`, split on the first `_`); a row whose half is empty is unannotated
  and gets dropped by the Hugo-symbol filter. `TotalFragments` no longer exists —
  `Tumor_Variant_Count` is `SplitFrags + RealignedFrags + DiscordantFrags`.
  `TranscriptUp/Down` and `ExonUp/Down` populate `Site1/2_Ensembl_Transcript_Id`
  and `Site1/2_Region_Number`. The converter is 3.0-only; 2.3 CSVs are rejected
- `Site1/2_Region` uses cBioPortal's vocabulary `{5_Prime_UTR, 3_Prime_UTR, Promoter,
Exon, Intron}`. Isofox records a transcript + exon rank only when the breakend falls
  on an exon (`FusionReadData.java` `TransExonRef`), so: transcript + exon → `Exon`;
  gene named but no exon ref → `Intron` (the caller assigns genes by gene-body overlap
  with no promoter allowance); no gene → `NA`. UTR/Promoter are not derivable —
  `pass_fusions.tsv` has no coding-type or CDS boundaries, and UTR bases are exonic
- `ml_format_cnv.R` / `ml_format_expression.R` check `basename(input)` — inputs must be named `data_cna_long.txt` / `data_expression.txt`
- No internet on compute nodes — `NXF_OFFLINE=true`; pre-pull containers on login nodes
- VEP/PCGR data must be pre-staged

## Timeline Generation (`gen_timeline.R`)

Generates a single combined `data_timeline.txt` from ARGO clinical CSVs. All event types are merged into one file, distinguished by EVENT_TYPE column:

- `SURGERY` — surgical treatments (subtype, site, intent)
- `TREATMENT` — systemic therapies (drug, type, intent)
- `STATUS` — follow-up disease status
- `SPECIMEN` — specimen collection (site, type, sample type)
- `LAB_TEST` — biomarker results (PSA, ER/PR/HER2, CEA, etc.)

Columns are the union of all event types; columns not applicable to a given EVENT_TYPE are empty. Common columns first (PATIENT_ID, START_DATE, STOP_DATE, EVENT_TYPE), then the rest alphabetically.

Key behaviors:

- **START_DATE defaults to 0 when missing** — any NA/empty date becomes day 0 (diagnosis day)
- Filters patients from `sample_registrations`: `Total DNA` + `Tumour` + `Solid tissue`
  rows matching `-\d+[A-Z]*[DR]T$`. `clin_format.R` applies the identical rule, so
  germline normals (buffy coat included) never become rows in `data_clinical_sample.txt`
  and a donor with no sequenced tumour drops out of every clinical output
- Every event of an eligible patient is kept, whichever specimen, treatment or diagnosis
  it hangs off — timeline events are only ever filtered by patient
- When `genomic_subjects` TSV is provided (both mode), restricts timeline output to genomic subjects only
- Nearest follow-up visit determines specimen SAMPLE_TYPE (Primary / Recurrence / Metastasis)

## Primary Diagnosis Selection (`clin_format.R`)

A donor may have several primary diagnoses. Each sample row takes the diagnosis of its own
sequenced specimen (`specimens.csv` → `submitter_primary_diagnosis_id`), so the germline
normal's diagnosis can never win; when the specimen carries no usable link the patient's
first diagnosis is used. Every other clinical table (treatments, surgeries, systemic
therapies, radiations, follow-ups, biomarkers) is merged by patient alone — patient-level
information is kept even when it is attached to a different diagnosis or specimen.

## MOHCCN Mapping Tables

Three mapping CSVs in `assets/` translate MOHCCN plain-text values to ontology/ICD-O codes:

- `mohccn_*_primary_site.csv` — ICD-O topography codes (reversed: code → label for surgery/specimen site)
- `mohccn_*_specimen_tissue_source.csv` — tissue source codes
- `mohccn_*_treatment_intent.csv` — treatment intent codes

Params: `mohccn_primary_site_map`, `mohccn_specimen_tissue_source_map`, `mohccn_treatment_intent_map`. Used by both `clin_format.R` and `gen_timeline.R`.

## Incremental Processing

The genomic workflow checks for pre-existing output files per subject. If all expected outputs (CNV, SV, expression, mutations) already exist, processing is skipped. This allows adding new subjects to the samplesheet and re-running without reprocessing the entire cohort. Delete a subject's output directory to force reprocessing.

## Combining Two Runs (`bin/combine_cbioportal_outputs.py`)

Standalone CLI utility (not called by the pipeline) that merges two study folders or
`.tar.gz` archives — used to add a batch of samples to an already-loaded study.

Merge strategy per file category:

- row-append — `data_cna_hg38.seg`, `data_cna_long.txt`, `data_sv.txt`, `util_linking_file.txt`
  (the linking file has no trailing newline, so a separator is inserted)
- 2-line header — `data_mutations_dna_rna_germline.txt` (`#version 2.4` + column header)
- wide matrix — `data_expression.txt`, the six `data_mutational_signatures_*` files
- union-append — `data_timeline.txt`; its column set varies with the event types a batch
  contained, so columns are unioned in `gen_timeline.R` order (4 common, then alphabetical)
- clinical — `data_clinical_{sample,patient}.txt`; `clin_format.R` drops columns whose source
  data is absent, so columns are unioned and the 4 `#` metadata rows rebuilt per column
- headerless dedupe — `cancer_type.txt`
- discovered by glob — `meta_*.txt` and `case_lists/*.txt`, so new files are never silently lost
- union-copied — per-subject folders, keeping the merged folder usable for incremental resume

`machine_learning/` is deliberately **not** merged (log2/standardised tables are
cohort-normalised) — regenerate it over the combined cohort. Anything else is reported as
unhandled; `--strict` makes that a non-zero exit. Add new outputs to a dispatch table in the
script and cover them in `tests/test_combine_cbioportal_outputs.py`.

Overlapping sample IDs are a hard error. Patient rows are appended without deduplication.

## Study Archive

`PACKAGE_CBIOPORTAL` is invoked from `main.nf` (not from `workflows/genomic.nf`) so that in
`both` mode `<study>.tar.gz` mirrors the loadable part of the study directory — genomic +
clinical + timeline + `util_linking_file.txt` + case lists. `GENOMIC` and `CLINICAL` each
emit `package_files`; `main.nf` mixes them and calls the module once. Clinical-only mode
produces no archive (it has no `study_id`/group directory).

The archive's md5 is not reproducible (tar embeds per-file mtimes), so `*.tar.gz` is listed
in `tests/.nftignore` — assert on its contents, not its checksum.

## Process Labels (`conf/base.config`)

Default: 1 CPU, 1 GB, 4 min (scaled by `task.attempt`). See `conf/base.config` for full table.

## Mutational Signatures

See `docs/mutational_signatures.md` for detailed SBS/DBS/ID documentation.

## Testing

nf-test gotchas:

- Use `path(f.toString())` — channel file outputs are `String`, not `Path`
- Sort snapshots: `.sort { it.toString().split('/').last() }`
- `collectFile` with `storeDir` won't create dirs — call `file("${params.outdir}/GROUP").mkdirs()` in test setup
- `genomic_ml` uses `options "-stub-run"` to skip `DOWNLOAD_KNOWN_FUSIONS`
