# CLAUDE.md — ampliseq

Nextflow DSL2 pipeline converting Ampliseq VCFs + TSV exports + clinical files → cBioPortal-ready files. HPC-only (SLURM + Apptainer — **no Docker, no sudo**).

---

## Layout

```
workflows/ampliseq-cbioportal.nf    # Full DAG: per-sample → merge → deanon → clinical → meta
modules/local/                       # One process per file (see Data Flow below)
bin/                                 # Python scripts; also runnable standalone via run_pipeline.sh
assets/samplesheet.csv               # Full test samplesheet (SAMPLE_001 + SAMPLE_002)
assets/samplesheet_sample1.csv       # First-run input for incremental test (SAMPLE_001 only)
```

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

---

## Key Implementation Notes

- CNA copy-number mapping: `0→-2, 1→-1, 3→1, ≥4→2`; CN=2 is normal and dropped.
- `data_cna.txt` is long format (`Hugo_Symbol, Sample_Id, Value`); `meta_cna.txt` uses `datatype: DISCRETE_LONG`.
- vcf2maf container mounts `vep_data` as `/home/jbellavance/` inside the container.
- All deanon scripts warn on unmatched IDs but leave them unchanged.
- `filter_tsv_variants=true` (default): mutations filtered to TSV coordinates. `false`: all MAF rows pass through.
- `vcf_to_seg.py` sets `num.mark=1` (ampliseq VCFs carry no probe-count).
- `meta_seg.txt`: `datatype: SEG`, `show_profile_in_analysis_tab: false`.

---

## Incremental Runs

A subject is skipped if all four per-sample files exist under `{outdir}/samples/{sample_id}/`:
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
