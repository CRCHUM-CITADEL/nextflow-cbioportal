# CLAUDE.md — oncoanalyser_3.0.0

Oncoanalyser WGS/WTS + clinical CSVs → cBioPortal. HPC-only (SLURM + Apptainer).

**Supports oncoanalyser 3.0 only** (since 4.0.0). The 2.3 output layout is not
accepted — converters hard-require the 3.0 columns and `stop()` otherwise.
Requires Nextflow >= 26.04.4.

## Key Rules

- R scripts → `container_r`; Python → `container_python`; SigProfiler → `container_sigprofiler`;
  VCF→MAF → `container_mafsmith`
- **VCF→MAF is `MAFSMITH` (nf-osi/mafsmith), not vcf2maf** (since 4.0.0).
  `mafsmith vcf2maf` always reads `$HOME/.mafsmith`, ignoring `--data-dir`, so the
  module sets `HOME=./` and symlinks `mafsmith_data` in as `.mafsmith`.
  `genome_reference` is still required, but now only by SigProfiler
- `mafsmith_data` is optional: empty means `DOWNLOAD_MAFSMITH` runs `mafsmith fetch`
  (`storeDir assets/mafsmith`), gated on `ch_has_work` like the VEP/PCGR downloads.
  **Prefer pre-staging it.** The fetch pulls Ensembl's primary assembly (contigs
  `1`/`2`/`X`) while oncoanalyser aligns against GATK hg38 (`chr1`/`chr2`/`chrX`), and
  mafsmith's reference must match the genome the mutations were called on. A hand-built
  bundle takes `--gff3` and `--ref-fasta` — **both or neither**, since `fetch` only
  takes the link path when it has the pair, and otherwise downloads both
- fastvep is **built into the container**, not fetched at run time: `resolve_fastvep()`
  falls back to `which fastvep` when `<data-dir>/bin/fastvep` is missing. That is why
  `DOWNLOAD_MAFSMITH` can pass `--skip-fastvep` and the image needs no cargo toolchain
  at run time. Rebuilding `containers/mafsmith_v0.1.0.def` is required for this
- The `MAFSMITH` script locates the VCF sample columns by the `FORMAT` header and the
  MAF `Tumor_Sample_Barcode` column by name, rather than at fixed positions
- mafsmith writes **53 MAF columns to vcf2maf's ~133**. `mafsmith_retain_ann` names VEP
  CSQ subfields to append as extra columns; names must match fastVEP's `DEFAULT_CSQ_FIELDS`
  (an unknown name silently yields an empty column). `RefSeq` and `VARIANT_CLASS` have no
  fastVEP counterpart and cannot be restored
- `CONVERT_CPSR_TO_MAF` merges CPSR germline calls (`Mutation_Status=Germline`, filtered to
  Pathogenic / Likely_Pathogenic / VUS on the CPSR TSV's 52nd field) into the somatic+RNA MAF,
  then drops `Intron`/`IGR` rows. `gen_convert_cpsr_to_maf.R` keys each germline row on the
  **incoming MAF's header** and pins it back to those columns before appending — `maf_entry$X <- v`
  appends a slot when `X` is absent, and `rbind()` onto a zero-row frame widens silently rather
  than erroring, which writes rows wider than the header. Add a column here only if the annotator
  actually emits it
- `data_sv.txt` rows require Hugo symbols at both sites — filter unannotated rows
- SV classification (`gen_esvee_sv_to_cbioportal.R`): BND ALT strand → `(+,-)` DEL, `(-,+)` DUP, `(+,+)/(−,−)` INV, diff chr TRANSLOC. DNA SVs: `DNA_Support=Yes, RNA_Support=No`
- RNA fusions (`gen_isofox_fusion_to_cbioportal.R`): `Class=FUSION, DNA_Support=No, RNA_Support=Yes`. Both merge into `data_sv.txt`
- Isofox 3.0 `pass_fusions.tsv` packs both gene symbols into one `Name` column
  (`<up>_<down>`, split on the first `_`); a row whose half is empty is unannotated
  and gets dropped by the Hugo-symbol filter. `TotalFragments` no longer exists —
  `Tumor_Variant_Count` is `SplitFrags + RealignedFrags + DiscordantFrags`.
  `TranscriptUp/Down` and `ExonUp/Down` populate `Site1/2_Ensembl_Transcript_Id`
  and `Site1/2_Region_Number`. The converter requires the 3.0 columns and rejects 2.3 input
- `Site1/2_Region` uses cBioPortal's vocabulary `{5_Prime_UTR, 3_Prime_UTR, Promoter,
Exon, Intron}`. Isofox records a transcript + exon rank only when the breakend falls
  on an exon (`FusionReadData.java` `TransExonRef`), so: transcript + exon → `Exon`;
  gene named but no exon ref → `Intron` (the caller assigns genes by gene-body overlap
  with no promoter allowance); no gene → `NA`. UTR/Promoter are not derivable —
  `pass_fusions.tsv` has no coding-type or CDS boundaries, and UTR bases are exonic
- `ml_format_cnv.R` / `ml_format_expression.R` check `basename(input)` — inputs must be named `data_cna_long.txt` / `data_expression.txt`
- No internet on compute nodes — `NXF_OFFLINE=true`; pre-pull containers on login nodes
- VEP/PCGR data must be pre-staged
- Nextflow optional outputs: `optional: true` is an option on the whole output
  declaration, NOT an argument to `path()`. `tuple val(x), path("f.txt"), emit: y,
optional: true` works; `path("f.txt", optional: true)` silently does nothing and
  the task fails with `MissingFileException` when the file is genuinely absent
- New `bin/` scripts must be `chmod +x`, or the process fails with exit 126
  ("Permission denied") inside the container rather than a normal error

## Configuration Layout

`nextflow.config` carries **portable defaults only** — public `oras://` container
images, empty reference-data paths, `outdir = 'output'`. It cannot be renamed:
Nextflow auto-loads `$projectDir/nextflow.config`, and it is the only home of the
`test` profile (which `nf-test.config` selects), the nf-schema plugin declaration,
the manifest, and the `conf/*.config` includes.

Site-specific values live in `nextflow_citadel.config`, loaded by the `citadel`
profile — `/project/60005` reference data, prebuilt `.sif` containers,
`--account=def-chasse`, scratch, and the shared Apptainer cache:

```bash
nextflow run main.nf -profile slurm,apptainer,citadel ...   # CRCHUM; citadel LAST
nextflow run main.nf -profile apptainer --ensembl_annotations ... --genome_reference ...
```

`ensembl_annotations` and `genome_reference` have no default and are required for
genomic/both mode; `mafsmith_data`, `vep_data` and `pcgr_data` are optional and
auto-download when empty. (`ensembl_annotations_expr` was dropped in 4.0.0 —
Isofox expression now maps Ensembl→Entrez through `ensembl_annotations` too.) Note
`.gitignore` has a `/nextflow_*.config` rule with an explicit negation for
`nextflow_citadel.config`.

`--incremental` marks a follow-up load into a study cBioPortal already holds: it
skips `cancer_type.txt` / `meta_cancer_type.txt` so the cancer type is not registered
twice. It is unrelated to per-subject output caching, which is automatic.

## Clinical Module Layout

Clinical output is built by small, focused modules so a failure names the output
it was building. `bin/clinical_common.R` holds the helpers both halves share
(`parse_day_interval`, `clean_json_array`, `read_mohccn_map`, `reverse_mohccn_map`,
`apply_mohccn_map`, `apply_ps_label`, `is_cbio_sample`, `read_genomic_subjects`,
`rbind_fill`) and is staged into each process as an explicit `path()` input — not
relied on via `$projectDir/bin` on PATH — so it participates in the task hash and
`-resume` invalidates correctly when it changes.

```
BUILD_CLINICAL_TABLE → clinical_merged.tsv
   ├→ WRITE_CLINICAL_SAMPLE  → data_clinical_sample.txt
   └→ WRITE_CLINICAL_PATIENT → data_clinical_patient.txt

GENERATE_TIMELINE_{SURGERY,TREATMENT,STATUS,SPECIMEN,LAB_TEST} → MERGE_TIMELINE
   → data_timeline.txt
```

- `BUILD_CLINICAL_TABLE` runs the loaders + merge cascade **once** per group and
  writes an inspectable intermediate table; the two writers only apply `col_defs`.
  (Previously one `FORMAT_CLINICAL` ran twice per group, repeating ~95% of the work
  just to emit a different column list.)
- The intermediate is written with `na="NA"` and read back with
  `colClasses="character"` so nothing is reformatted on the round trip — that is
  what keeps the split byte-identical.
- Each per-event timeline script writes nothing when it has no rows; its output is
  optional, so `MERGE_TIMELINE` routinely sees a subset of the five parts.
- `MERGE_TIMELINE` owns `rbind_fill` + the column order (4 common columns, then
  alphabetical) — a hard contract with `combine_cbioportal_outputs.py`. It also
  stable-sorts rows by a fixed EVENT_TYPE order (SURGERY, TREATMENT, STATUS,
  SPECIMEN, LAB_TEST), because `mix()`/`groupTuple()` guarantees no arrival order
  and the output would otherwise be nondeterministic.

## Timeline Generation (`gen_timeline_*.R` → `merge_timeline.R`)

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
  rows matching `-\d+[A-Z]*[DR]T$`. The shared `is_cbio_sample()` in `clinical_common.R`
  is the single implementation of this rule, so
  germline normals (buffy coat included) never become rows in `data_clinical_sample.txt`
  and a donor with no sequenced tumour drops out of every clinical output
- Every event of an eligible patient is kept, whichever specimen, treatment or diagnosis
  it hangs off — timeline events are only ever filtered by patient
- When `genomic_subjects` TSV is provided (both mode), restricts timeline output to genomic subjects only
- Nearest follow-up visit determines specimen SAMPLE_TYPE (Primary / Recurrence / Metastasis)

## Primary Diagnosis Selection (`build_clinical_table.R`)

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

Params: `mohccn_primary_site_map`, `mohccn_specimen_tissue_source_map`, `mohccn_treatment_intent_map`. Used by `build_clinical_table.R` and the `gen_timeline_*.R` scripts, all through
`clinical_common.R`'s `read_mohccn_map()` / `reverse_mohccn_map()`.

Three sample-level columns carry an ICD-O code in the source data and get a
reverse-mapped `*_LABEL` sibling via `apply_ps_label()`: `RELAPSE_SITE`,
`CANCER_TYPE_CODE` and `TUMOR_TISSUE_SITE`. Note `apply_ps_label()` resolves to the
**C+2-digit parent** code, because the MOHCCN primary-site table only carries parent
codes — `C22.0` and `C22.1` both label as "Liver and intrahepatic bile ducts".
Sub-site precision would need a fuller ICD-O-3 topography table as a new asset.
The `_LABEL` suffix is deliberate: cBioPortal reserves `CANCER_TYPE` /
`CANCER_TYPE_DETAILED` for OncoTree codes and would misread ICD-O labels there.

## Incremental Processing

Re-running with a samplesheet of old + new subjects skips the subjects already
processed and merges their published outputs with the freshly computed ones.

- **`subjectOutputs()` in `workflows/genomic.nf` is the single source of truth** for
  what a run publishes per subject. Both the "is this subject done" gate and the cache
  loader read that one table, and `cachedFiles()` asserts its tag argument is in it.
  Add a new per-subject output there and nowhere else — keeping the gate and the loader
  as two hand-written lists is how the ID signatures ended up gated but never reloaded,
  silently dropping every cached subject from the two `*_ID.txt` files.
- **A subject is cached all-or-nothing.** `resolveSubjectCache()` returns the whole set
  or `null`; a subject missing any `required: true` output is reprocessed from source
  and contributes _nothing_ from disk. Mixing stale and fresh files for one sample is
  not just duplicate rows — fresh and cached files share a basename, so the staged-list
  merges (`MERGE_EXPRESSION_FILES_TO_CBIOPORTAL`, `MERGE_SIGS_*`) abort the run with an
  input file name collision.
- The signature outputs are `required: false`: SigProfiler legitimately emits nothing
  for a subject with too few mutations, so gating on them would reprocess that subject
  forever. They are still always loaded when present.
- **Cached subjects bypass `MERGE_SAMPLE_SV`.** The published `<sample>.data_sv.txt` is
  already that module's merged DNA + RNA-fusion output, and re-merging it would hand the
  process an input with the same name as its own output.
  `<sample>.isofox_fusion.data_sv.txt` is never published, so there is nothing to cache.
- **The samplesheet defines the cohort.** The group-level files are rebuilt from it on
  every run, so a subject dropped from the samplesheet vanishes from the study even
  though its folder stays on disk. The workflow warns about this by diffing the previous
  run's `util_linking_file.txt`, read before it is overwritten.
- Delete a subject's output directory to force reprocessing.
- With no subject left to process, no reference data is downloaded: `GENOMIC_MUTATIONS`
  gates `DOWNLOAD_VEP_TEST` / `DOWNLOAD_PCGR` on `som_dna_vcf.count().filter { it > 0 }`.
- Case-list sample lists use `toList()`, not `collect()`: `collect()` emits nothing for
  an empty channel, which would shift the positional zip in `GENERATE_CASE_LIST` and
  write one modality's samples under another's label.

## Combining Two Runs (`bin/combine_cbioportal_outputs.py`)

Standalone CLI utility (not called by the pipeline) that merges two study folders or
`.tar.gz` archives — used to add a batch of samples to an already-loaded study.

Merge strategy per file category:

- row-append — `data_cna_hg38.seg`, `data_cna_long.txt`, `data_sv.txt`, `util_linking_file.txt`
  (the linking file has no trailing newline, so a separator is inserted)
- 2-line header — `data_mutations_dna_rna_germline.txt` (`#version 2.4` + column header)
- wide matrix — `data_expression.txt`, the six `data_mutational_signatures_*` files
- union-append — `data_timeline.txt`; its column set varies with the event types a batch
  contained, so columns are unioned in `merge_timeline.R` order (4 common, then alphabetical)
- clinical — `data_clinical_{sample,patient}.txt`; `write_clinical_table.R` drops columns whose source
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
- A process with no `stub:` block runs its **real** script under `-stub-run` — Nextflow
  does not error. Every module on the genomic path therefore needs one
- `tests/incremental_pipeline.nf.test` drives the whole pipeline under `-stub-run`:
  freshly processed subjects emit empty placeholders while cached subjects contribute
  real fixture content, so a sample's rows appear in a merged file if and only if it was
  read from the cache. Its config disables the container engines, since a stub is plain
  shell and pulling an image for a `touch` only invites registry-auth failures
- `WorkflowTask` exposes the public fields `name` and `success` only — use
  `workflow.trace.succeeded().collect { it.name }`; it has no `toString()`
