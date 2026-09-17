---
name: oncoanalyser migration
description: Pipeline was refactored from DRAGEN inputs to nf-core/oncoanalyser outputs; summarises what changed and design decisions
type: project
---

Pipeline was refactored to accept nf-core/oncoanalyser outputs instead of DRAGEN outputs.

**Why:** User requested oncoanalyser compatibility so the pipeline can consume HMFtools outputs (PAVE, SAGE, PURPLE, ESVEE, Isofox).

**How to apply:** All new development should assume the oncoanalyser **3.0** output directory structure. As of 4.0.0 the 2.3 layout is not supported — the converters hard-require the 3.0 columns and `stop()` otherwise.

## Source mapping

Paths below are relative to the samplesheet's `folder` column, and `<subject>` is `subject_id`. These are what `workflows/genomic.nf` actually resolves — keep this table in sync with the `findOncoFile` calls there.

| Modality         | Old (DRAGEN)                                      | New (oncoanalyser 3.0)                                   |
| ---------------- | ------------------------------------------------- | -------------------------------------------------------- |
| Mutations        | `*.WGS_somatic-tumor_normal.hard-filtered.vcf.gz` | `pave/somatic/<subject>-T.pave.somatic.vcf.gz`           |
| Germline         | `*.WGS_germinal.hard-filtered.vcf.gz`             | `pave/germline/<subject>-T.pave.germline.vcf.gz`         |
| RNA append       | N/A                                               | `sage_append/<subject>-T/<subject>-T.sage.append.vcf.gz` |
| CNV (segments)   | `*.WGS_somatic-tumor_normal.cnv.vcf.gz`           | `purple/<subject>-T.purple.cnv.somatic.tsv`              |
| CNV (genes)      | derived from VCF + annotation                     | `purple/<subject>-T.purple.cnv.gene.tsv`                 |
| SV (DNA)         | N/A                                               | `esvee/<subject>-T.esvee.somatic.vcf.gz`                 |
| Fusions (RNA)    | `*.fusion_candidates.final`                       | `isofox/<subject>-T.isf.pass_fusions.tsv`                |
| Expression       | `*.quant.genes.sf` (Salmon)                       | `isofox/<subject>-T.isf.gene_data.tsv`                   |
| Signature counts | N/A                                               | `sigs/<subject>-T.sig.snv_counts.csv`                    |

## Samplesheet

One row per subject (`assets/schema_genomic_input.json`):
`group, subject_id, sample_id, folder`

`folder` is the per-subject oncoanalyser output directory; every modality file above is resolved relative to it. A missing file logs a warning and that modality is skipped for the subject, rather than failing the run.

## Still present (do not assume removed)

An earlier version of this note claimed PCGR/CPSR and `container_pcgr` had been dropped. They have not:

- `modules/local/pcgr/` and `modules/local/convert_cpsr_to_maf/` both exist and are used
- `container_pcgr` is a live param
- RNA variant integration into the DNA MAF also still exists (`modules/local/integrate_rna_variants/`)

## New files

- `bin/gen_purple_cnv_to_cbioportal.R` — PURPLE TSVs → SEG + DISCRETE_LONG
- `bin/gen_esvee_sv_to_cbioportal.R` — ESVEE BND VCF → cBioPortal data_sv.txt (gene annotation via foverlaps)
- `bin/gen_isofox_fusion_to_cbioportal.R` — Isofox `pass_fusions.tsv` → cBioPortal data_sv.txt
- `bin/gen_isofox_expression_to_cbioportal.R` — Isofox `gene_data.tsv` → cBioPortal expression format
- `modules/local/purple_cnv_to_cbioportal/`, `esvee_sv_to_cbioportal/`, `isofox_fusion_to_cbioportal/`, `isofox_expression_to_cbioportal/`
