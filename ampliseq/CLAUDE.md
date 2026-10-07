# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Nextflow DSL2 pipeline converting Ampliseq VCFs + TSV exports + clinical files → cBioPortal-ready files. HPC-only (SLURM + Apptainer — **no Docker, no sudo**).

---

## Commands

```bash
# Run tests (from ampliseq/ directory) — local ./nf-test binary; `+` appends to nf-test.config's `test` profile
./nf-test test tests/modules --profile=+apptainer        # all module tests; MAFSMITH pulls the private ghcr image
./nf-test test tests/modules/package_cbioportal.nf.test --profile test   # what CI runs (no containers)

# Run locally with test data (skips mafsmith)
nextflow run main.nf -profile test,apptainer

# Run with real data
nextflow run main.nf -profile apptainer \
  --input samplesheet.csv \
  --outdir results/ \
  --patient_file patient_file.txt \
  --sample_file sample_file.txt \
  --linking_file linking_file.txt \
  --mafsmith_data /path/to/mafsmith_home/ \
  --study_id my_study

# Resume a failed/interrupted run
nextflow run main.nf ... -resume
```

The `test` profile sets `skip_vcf2maf=true` and uses stub MAFs from `assets/`, avoiding the need for the mafsmith reference bundle.

Nextflow: **26.04.4 minimum** (`nextflowVersion '!>=26.04.4'`); the ampliseq CI job overrides `NXF_VER` to 26.04.4 (the workflow-wide default is still 25.10.2 for the other pipelines). See **Nextflow 26 notes** below before touching params or `publishDir`.

---

## Workflow Architecture

The main DAG in `workflows/ampliseq-cbioportal.nf` orchestrates 3 subworkflows and 2 standalone modules:

```
ch_samplesheet
  │
  ├─ branch → existing (skip) / new_sample
  │
  ├─ PER_SAMPLE_FORMAT (subworkflow) ── new samples only
  │    ├─ FORMAT_SV         (TSV → _sv.txt)
  │    ├─ FORMAT_CNA        (TSV → _cna.txt)
  │    ├─ VCF_TO_SEG        (CNV VCF → _seg.txt)
  │    └─ MAFSMITH (+ DOWNLOAD_MAFSMITH) or STUB_MAF → FILTER_MUTATIONS or PASSTHROUGH_MUTATIONS → _mutations.txt
  │
  ├─ FILTER_LINKING (module) ── restrict linking file to samplesheet samples
  │
  ├─ MERGE_DEANON (subworkflow) ── merge all per-sample files + deanonymise
  │    ├─ MERGE_{SV,CNA,MUTATIONS,SEG}
  │    └─ DEANON_{MUTATIONS,SV,CNA,SEG}
  │
  ├─ STUDY_METADATA (subworkflow) ── clinical files + meta + case lists
  │    ├─ CLINICAL_PATIENTS / CLINICAL_SAMPLES
  │    └─ WRITE_CASE_LISTS / WRITE_META
  │
  └─ PACKAGE_CBIOPORTAL (module) ── tar.gz all outputs for transfer
```

**Key branching logic:** `params.skip_vcf2maf` chooses between real VCF→MAF conversion (mafsmith + fastVEP, GRCh37, PASS records only) and stub MAFs. `params.filter_tsv_variants` (default false) chooses between FILTER_MUTATIONS, which keeps MAF rows whose `[Start_Position, End_Position]` overlaps any `analysis_*_export.tsv` row's `Chr:Start-End` (inclusive, `chr` prefix ignored, all variant types), and passing all MAF rows through. It used to require an exact `Start-End` match, which dropped every SNV because TSV rows are regions.

---

## Input Files (per-sample folder)

- `analysis_*_export.tsv` — columns: `Chr, Start, End, Variant Type, Variant Subtype, Genes, Breakend Genes, Supporting Reads, Copy Number`
- `*-basespace-pisces.final.vcf.gz` — somatic mutations VCF; filename prefix = `SAMPLE_ID`
- `*-basespace-cnv.final.vcf` — CNV VCF; needs `CN` in FORMAT and `END` in INFO

**Patient file** columns: `patient_id, moh_id, age, sex, os_status, os_months, smoking_history` (`moh_id` → `MOH_ID` in `data_clinical_patient.txt`; must be kept).
**Sample file** columns: `sample_id, patient_id, cancer_type, cancer_type_detailed, sample_type, primary_tumor_site, metastatic_tumor_site, tumor_purity`.
`bin/clinical_{patients,sample}_format.py` index these by name, so a missing column is a `KeyError`. Extra columns are ignored. Keep `assets/{patient,sample}_file.txt` in sync with the scripts.

**Linking file** (`linking_file.txt`): maps anonymized → real IDs.
```
sample_id   deanon_sample_id   deanon_patient_id
```
The `sample_id` column in the **sample file** must use deanonymized IDs (`deanon_sample_id`), not anonymized ones.

---

## Data Flow

**Per-sample** → published to `{outdir}/samples/{sample_id}/`:
1. `analysis_*_export.tsv` → FORMAT_SV → `_sv.txt` (`Variant Subtype = FUSION`)
2. `analysis_*_export.tsv` → FORMAT_CNA → `_cna.txt` (`DUPLICATION`/`DELETION`)
3. VCF → MAFSMITH (mafsmith + fastVEP, GRCh37, PASS only; header on line 1, `chr` prefix kept) → FILTER_MUTATIONS or PASSTHROUGH_MUTATIONS → `_mutations.txt`
4. `*-cnv.final.vcf` → VCF_TO_SEG → `_seg.txt` (PASS only; `seg.mean = log2(CN/2)`; CN=0 → −3.0)

**Downstream** (re-runs on every execution over all samples):
5. FILTER_LINKING → linking filtered to samplesheet samples only
6. MERGE → DEANON → `data_mutations.txt`, `data_sv.txt`, `data_cna.txt`, `data_seg.txt`
7. CLINICAL_PATIENTS / CLINICAL_SAMPLES → filtered to samplesheet patients/samples
8. WRITE_CASE_LISTS + WRITE_META
9. PACKAGE_CBIOPORTAL → `{study_id}.tar.gz`

---

## Key Implementation Notes

- CNA copy-number mapping: `0→-2, 1→-1, 3→1, ≥4→2`; CN=2 is normal and dropped.
- `data_cna.txt` is long format (`Hugo_Symbol, Sample_Id, Value`); `meta_cna.txt` uses `datatype: DISCRETE_LONG`.
- All deanon scripts warn on unmatched IDs but leave them unchanged.
- `vcf_to_seg.py` sets `num.mark=1` (ampliseq VCFs carry no probe-count).
- `meta_seg.txt`: `datatype: SEG`, `show_profile_in_analysis_tab: false`.

## Container Labels

Two process labels control container assignment in `nextflow.config`:
- `python` → `params.python_sif` (local Apptainer image built from `containers/python-ampliseq.def`)
- `mafsmith` → `params.mafsmith_container` (mafsmith + fastVEP, built from `containers/mafsmith-fastvep_v0.1.0-0.4.0.def`)

`params.mafsmith_data` (optional) is the mafsmith home with `GRCh37/{reference.fa,genes.gff3.gz}`. If it is null, empty or the string `"null"`, `DOWNLOAD_MAFSMITH` fetches it once into `assets/mafsmith` (storeDir). `MAFSMITH` fails early if either reference file is missing or empty, unless `ext.args` contains `--skip-annotation` (the module test relies on this). **Gotcha:** `-stub-run` with `skip_vcf2maf=false` stores an *empty* `assets/mafsmith`, which storeDir then reuses; delete it after stub runs. **Contig naming:** the bundle is assumed to use `mafsmith fetch` (Ensembl) naming: `1, 2, X, MT`. mafsmith refuses a VCF whose contigs are not in the FASTA (`VCF chromosome 'chr1' was not found in the FASTA index`); the GFF's naming does not matter (verified on real Ensembl GRCh37 chr12). `MAFSMITH` therefore strips `chr` from the PASS VCF (records + `##contig`, `chrM`→`MT`) and puts it back on the MAF `Chromosome` column (`MT`→`chrM`) when the source VCF was chr-prefixed, as Pisces VCFs are. A pre-staged `--mafsmith_data` with chr-named FASTA is not supported. `Entrez_Gene_Id` is always 0 with mafsmith/fastVEP (no Entrez source); only `Hugo_Symbol = Unknown` signals missing annotation. Mutation_Status stays blank (tumor-only). `MERGE_MUTATIONS` merges per-sample files by column name, so older vcf2maf-era files (different columns) still merge correctly.

Container definitions live in `containers/`. See `containers/README.md` for local build, `--mafsmith_container <local.sif>`, and the push/tag convention. ghcr images are private, so pulls need `apptainer registry login`.

---

## Nextflow 26 notes

- **CLI params are Strings on NF 26**: `--skip_vcf2maf false` arrives as `"false"`, which is truthy. Boolean params are read as `params.x.toString().toBoolean()` (`skip_vcf2maf`, `filter_tsv_variants`, `anonymize`). Do the same for any new boolean param.
- **`publishDir` paths that use `meta` must be closures**: `publishDir { "${params.outdir}/samples/${meta.sample_id}" }`. A plain GString fails on 26 with `No such variable: meta`.
- **nf-schema 2.8.0** (requires NF >= 26.04): 2.5.1 warns `Unrecognized config option 'validation.*'` on 26. From 2.7 on, path params in `nextflow_schema.json`/config must default to `null`, not `""`.

---

## Incremental Runs

A sample is skipped if all four per-sample files exist under `{outdir}/samples/{sample_id}/`:
`{id}_sv.txt`, `{id}_cna.txt`, `{id}_seg.txt`, `{id}_mutations.txt`

Merge/deanon/clinical steps always re-run over all samples combined. Use the same `--outdir` across runs.

---

## Standalone Scripts (`bin/`)

All Python scripts write output relative to `os.getcwd()` — run from the target output directory:
```bash
cd /path/to/output
python3 /path/to/bin/format_tsv.py    <export.tsv>  <SAMPLE_ID>
python3 /path/to/bin/format_cna.py    <export.tsv>  <SAMPLE_ID>
python3 /path/to/bin/vcf_to_seg.py    <cnv.vcf>     <SAMPLE_ID>
python3 /path/to/bin/format_mutations.py data_mutations.txt <linking_file>
python3 /path/to/bin/format_sv.py       data_sv.txt         <linking_file>
python3 /path/to/bin/format_cna_deanon.py data_cna.txt      <linking_file>
python3 /path/to/bin/seg_deanon.py     data_seg.txt         <linking_file>
```
`bin/run_pipeline.sh` orchestrates all of the above — contains **hardcoded cluster paths** that must be updated before use.
