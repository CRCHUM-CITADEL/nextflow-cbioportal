---
name: oncoanalyser migration
description: Pipeline was refactored from DRAGEN inputs to nf-core/oncoanalyser outputs; summarises what changed and design decisions
type: project
---

Pipeline was refactored to accept nf-core/oncoanalyser outputs instead of DRAGEN outputs.

**Why:** User requested oncoanalyser compatibility so the pipeline can consume HMFtools outputs (SAGE, PURPLE, ESVEE, Isofox).

**How to apply:** All new development should assume oncoanalyser output directory structure.

## Source mapping

| Modality | Old (DRAGEN) | New (oncoanalyser) |
|---|---|---|
| Mutations | `*.WGS_somatic-tumor_normal.hard-filtered.vcf.gz` | `sage/<tumor_id>.sage.vcf.gz` |
| CNV (segments) | `*.WGS_somatic-tumor_normal.cnv.vcf.gz` | `purple/<tumor_id>.purple.cnv.somatic.tsv` |
| CNV (genes) | derived from VCF + annotation | `purple/<tumor_id>.purple.cnv.gene.tsv` |
| SV (DNA) | N/A | `esvee/<tumor_id>.esvee.somatic.vcf.gz` |
| Fusions (RNA) | `*.fusion_candidates.final` | `isofox/<rna_id>.isofox.fusion.tsv` |
| Expression | `*.quant.genes.sf` (Salmon) | `isofox/<rna_id>.isofox.exp.tsv` |
| Germline | `*.WGS_germinal.hard-filtered.vcf.gz` → PCGR/CPSR | **Removed** |

## Samplesheet

Now identical to the oncoanalyser INPUT samplesheet:
`group_id, subject_id, sample_id, sample_type (tumor/normal), sequence_type (dna/rna), filetype, filepath`

Multiple rows per sample (one per file type); deduplicated in the workflow.

`oncoanalyser_outdir` param points to where oncoanalyser wrote its results.

## Removed components

- PCGR/CPSR germline annotation (no direct equivalent in oncoanalyser)
- RNA variant integration into DNA MAF (oncoanalyser uses separate Isofox for RNA)
- `container_pcgr` parameter

## New files

- `bin/gen_purple_cnv_to_cbioportal.R` — PURPLE TSVs → SEG + DISCRETE_LONG
- `bin/gen_esvee_sv_to_cbioportal.R` — ESVEE BND VCF → cBioPortal data_sv.txt (gene annotation via foverlaps)
- `bin/gen_isofox_fusion_to_cbioportal.R` — Isofox fusion.tsv → cBioPortal data_sv.txt
- `bin/gen_isofox_expression_to_cbioportal.R` — Isofox exp.tsv → cBioPortal expression format
- `modules/local/purple_cnv_to_cbioportal/`, `esvee_sv_to_cbioportal/`, `isofox_fusion_to_cbioportal/`, `isofox_expression_to_cbioportal/`
