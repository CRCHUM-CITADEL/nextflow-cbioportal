# CLAUDE.md — ampliseq

Ampliseq VCFs + TSV exports + clinical files → cBioPortal. HPC-only (SLURM + Apptainer).

## Structural Differences from oncoanalyser/dragen

- No `conf/` directory (config inline in `nextflow.config`), no `modules/nf-core/`
- Test data in `assets/` (no `test_data/` subdirectory)
- Container labels: `python` → `params.python_sif`; `vcf2maf` → `params.vcf2maf_container`

## Key Rules

- CNA copy-number mapping: `0→-2, 1→-1, 3→1, ≥4→2`; CN=2 is normal and dropped
- `data_cna.txt` is long format; `meta_cna.txt` uses `datatype: DISCRETE_LONG`
- `vcf_to_seg.py`: `num.mark=1`, `seg.mean` = raw CN (not log2); CN=2 rows dropped
- `params.skip_vcf2maf` chooses between real VCF→MAF (VEP v113, GRCh37/hg19) and stub MAFs
- `params.anonymize` enables clinical anonymization via BUILD_ANON_LINKING + ANONYMIZE_CLINICAL
- All deanon scripts warn on unmatched IDs but leave them unchanged
- Incremental: sample skipped if all 4 per-sample files exist in `{outdir}/samples/{sample_id}/`

## Workflow DAG

`PER_SAMPLE_FORMAT` (new samples) → `MERGE_DEANON` → `STUDY_METADATA` → `PACKAGE_CBIOPORTAL`

## Input Files (per-sample folder)

- `analysis_*_export.tsv` — Variant/CNA data
- `*-basespace-pisces.final.vcf.gz` — somatic mutations VCF
- `*-basespace-cnv.final.vcf` — CNV VCF (needs `CN` in FORMAT, `END` in INFO)
- `*-star-fusion.final.vcf` — (optional) fusion VCF

## Standalone Scripts

All `bin/` Python scripts write to `os.getcwd()` — run from target output dir.
`bin/run_pipeline.sh` has hardcoded cluster paths that must be updated.
