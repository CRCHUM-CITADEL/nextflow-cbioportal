# CLAUDE.md — dragen

DRAGEN somatic/germline + RNA + clinical CSVs → cBioPortal. Sister pipeline to `oncoanalyser/` — same structure and conventions. HPC-only (SLURM + Apptainer).

## Key Differences from oncoanalyser

- `data_sv.txt`: dragen fusion script (`gen_format_dragen_fusion.R`) does **not** filter unannotated rows (unlike oncoanalyser which does)
- Process labels have higher defaults — see `conf/base.config`
- No SigProfiler / mutational signatures

## Key Rules

- R scripts → `container_r`; Python → `container_python`
- `ml_format_cnv.R` / `ml_format_expression.R` check `basename(input)` — inputs must be named `data_cna_long.txt` / `data_expression.txt`
- No internet on compute nodes — `NXF_OFFLINE=true`; pre-pull containers on login nodes
- VEP/PCGR data must be pre-staged

## Testing

Same nf-test gotchas as oncoanalyser (see `oncoanalyser/CLAUDE.md`).
