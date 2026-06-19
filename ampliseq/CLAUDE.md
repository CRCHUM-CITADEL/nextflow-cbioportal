# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

Nextflow DSL2 pipeline converting Ampliseq VCFs + TSV exports + clinical files → cBioPortal-ready files. HPC-only (SLURM + Apptainer — **no Docker, no sudo**).

---

## Commands

```bash
# Run tests (from ampliseq/ directory)
nf-test test tests/modules/package_cbioportal.nf.test --profile test

# Run locally with test data (skips vcf2maf)
nextflow run main.nf -profile test,apptainer

# Run with real data
nextflow run main.nf -profile apptainer \
  --input samplesheet.csv \
  --outdir results/ \
  --patient_file patient_file.txt \
  --sample_file sample_file.txt \
  --linking_file linking_file.txt \
  --vcf2maf_container community.wave.seqera.io/library/vcf2maf_ensembl-vep:1b486a30e76e2908 \
  --vep_data /path/to/vep_data/ \
  --ref_fasta /path/to/hg19.fa \
  --study_id my_study

# Resume a failed/interrupted run
nextflow run main.nf ... -resume
```

The `test` profile sets `skip_vcf2maf=true` and uses stub MAFs from `assets/`, avoiding the need for VEP data and reference FASTA.

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
  │    └─ VCF_TO_MAF or STUB_MAF → FILTER_MUTATIONS or PASSTHROUGH_MUTATIONS → _mutations.txt
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

**Key branching logic:** `params.skip_vcf2maf` chooses between real VCF→MAF conversion (VEP v113, GRCh37/hg19) and stub MAFs. `params.filter_tsv_variants` chooses between filtering mutations to TSV coordinates or passing all MAF rows through.

---

## Input Files (per-sample folder)

- `analysis_*_export.tsv` — columns: `Chr, Start, End, Variant Type, Variant Subtype, Genes, Breakend Genes, Supporting Reads, Copy Number`
- `*-basespace-pisces.final.vcf.gz` — somatic mutations VCF; filename prefix = `SAMPLE_ID`
- `*-basespace-cnv.final.vcf` — CNV VCF; needs `CN` in FORMAT and `END` in INFO

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
3. VCF → VCF_TO_MAF (vcf2maf, VEP v113, GRCh37/hg19) → FILTER_MUTATIONS or PASSTHROUGH_MUTATIONS → `_mutations.txt`
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
- `vcf2maf` → `params.vcf2maf_container` (Wave-built image with vcf2maf + ensembl-vep)

The vcf2maf container mounts `vep_data` as `/home/jbellavance/` inside the container.

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
