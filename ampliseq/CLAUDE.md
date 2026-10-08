# CLAUDE.md — ampliseq

Ampliseq VCFs + TSV exports + clinical files → cBioPortal. HPC-only (SLURM + Apptainer).

## Structural Differences from oncoanalyser/dragen

- No `conf/` directory (config inline in `nextflow.config`), no `modules/nf-core/`
- Test data in `assets/` (no `test_data/` subdirectory)
- Container labels: `python` → `params.python_sif`; `vcf2maf` → `params.vcf2maf_container`

## Key Rules

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
  │    └─ MAFSMITH (+ DOWNLOAD_MAFSMITH → PREPARE_MAFSMITH) or STUB_MAF → FILTER_MUTATIONS or PASSTHROUGH_MUTATIONS → _mutations.txt
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

**Key branching logic:** `params.skip_vcf2maf` chooses between real VCF→MAF conversion (mafsmith + fastVEP, GRCh37, PASS records only) and stub MAFs. `params.filter_tsv_variants` (default false) turns on **filter mode**, which applies QC thresholds (all params, defaults shown) to three data types; with it off everything is kept:

- **Mutations** → FILTER_MUTATIONS instead of PASSTHROUGH_MUTATIONS: keep a MAF row only if its `[Start_Position, End_Position]` overlaps (inclusive, `chr` ignored) a TSV row with numeric `Depth >= mutation_min_depth` (250) and `VAF > mutation_min_vaf` (0.03). Depth/VAF come from the TSV. CNV/fusion rows have `N/A` there, so they are never match regions: a fusion row can span tens of Mb. Alleles are not compared.
- **CNA** → `format_cna.py --min-copy-number --confidence`: keep rows with `Confidence == cna_confidence` (HIGH, case-insensitive) and raw `Copy Number >= cna_min_copy_number` (6, before the half-up rounding). Every deletion is dropped and every kept gene is Value 2.
- **Fusions** → FORMAT_SV ignores any `*-star-fusion.final.vcf` and runs `format_tsv_to_sv.py --min-supporting-reads`: `Supporting Reads` summed per `(Genes, Breakend Genes)` pair (`N/A` = 0), and every row of a pair whose total is `>= sv_min_supporting_reads` (1000) is kept.
- A missing `Depth`/`VAF`/`Confidence` column fails the task in filter mode. `data_seg.txt` is never filtered.

---

## Input Files (per-sample folder)

- `analysis_*_export.tsv` — real export columns: `Chr, Start, End, Length, Variant Type, Variant Subtype, Ref, Alt, Genes, Cdot, Pdot, Exons, VAF, Confidence, Depth, Depth_Ref, Depth_Alt, Region, Effect, Germline Classification, Somatic Classification, Clinical Significance, In Report, Aggregated Normal Frequency, Copy Number, Included in TMB calculation, Possibly germline variant, Supporting Reads, Breakend Genes, Breakend Exon`. Missing values are the literal `N/A`; SNP rows are points (`Start == End`), CNV/fusion rows are regions. Small variants are `Variant Type = SNP` (subtype `SNP`/`INDEL`); fusions are `CNV`/`SV` + `FUSION`, one row per breakpoint. Scripts read columns by name; `assets/samples/*` uses this layout.
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

**Downstream** (re-runs on every execution over all samples): 5. FILTER_LINKING → linking filtered to samplesheet samples only 6. MERGE → DEANON → `data_mutations.txt`, `data_sv.txt`, `data_cna.txt`, `data_seg.txt` 7. CLINICAL_PATIENTS / CLINICAL_SAMPLES → filtered to samplesheet patients/samples 8. WRITE_CASE_LISTS + WRITE_META 9. PACKAGE_CBIOPORTAL → `{study_id}.tar.gz`

---

## Key Implementation Notes

- **SV columns:** both writers (`format_tsv_to_sv.py` for the export TSV, `fusion_vcf_to_sv.py` for `*-star-fusion.final.vcf`) emit exactly `cbio_sv.SV_COLUMNS`, the full cBioPortal SV layout modelled on oncoanalyser's writers (5'/3' sites, `Length` = `NA` for translocations, split/discordant counts, RNA/DNA support, `Event_Info`, `Annotation`, Nirvana `ANNOTATION` → `External_Annotation`). Add a column in `cbio_sv.py`, never in one writer only. STAR-Fusion `_1`/`_2` breakend records are paired into one row (each site's `EXON_NUM` → `Site*_Region_Number`). TSV fusion rows have one `Chr`, but `End` is the partner's position on _its_ chromosome (ROS1 chr6 → CD74 `149784294` is chr5): `format_tsv_to_sv.py --gene-loci` places it with `assets/grch37_gene_loci.tsv.gz` (`Name`/chrom/start/end of every `gene` record in the mafsmith GRCh37 `genes.gff3.gz`), taking the `Breakend Genes` locus containing `End` (±10 kb), else the gene's only chromosome, else `NA` + warning. `Class` is always `FUSION`, `Exons`/`Breakend Exon` → `Site1/2_Region_Number`.
- **Merges** (`MERGE_{MUTATIONS,SV,CNA,SEG}`) all use `bin/merge_tsv_by_header.sh`, which maps rows onto the union of columns **by name** and tolerates header-only and 0-byte files. Never concatenate positionally: per-sample files from older runs or different writers have different columns. The old inline awk dropped a whole sample's rows when any input file was empty.
- **Mutation enrichment:** fastVEP has no dbSNP/ClinVar/COSMIC/1000G, so `MAFSMITH` runs `bin/maf_add_vcf_annotations.py`, which fills `dbSNP_RS`, `Existing_variation`, `CLIN_SIG` and `AF` from the VCF ID column and Nirvana `clinvar`/`cosmic`/`AF1000G` (rows matched on chrom + MAF-style start + alt; indels included). Tumor-only: `Matched_Norm_Sample_Barcode` and `Match_Norm_Seq_Allele1/2` are blanked (mafsmith writes `NORMAL` + the ref allele).
- **Deanon scripts** read with `dtype=str, keep_default_na=False`, so values (including `NA`) round-trip unchanged and numeric-looking IDs don't crash `.str.upper()`.
- **Case lists** carry `case_list_category`; `cases_cnaseq.txt` (`all_cases_with_mutation_and_cna_data`) is written for OncoPrint. `PACKAGE_CBIOPORTAL` copies `case_lists` with `cp -rL`: it is staged as a symlink, and plain `cp -r` packed a dangling link instead of the files.
- CNA copy-number mapping: the copy number is rounded half-up, then `≤0→-2, 1→-1, 3→1, ≥4→2`; CN=2 is normal and dropped. Multi-gene `Genes` cells (`,`/`;`) give one row per gene; duplicate genes keep the most extreme value.
- `vcf_to_seg.py` reads `CN` by its FORMAT key (FORMAT may be `CN` or e.g. `GT:CN`) and accepts `.vcf.gz`. **Open question:** it writes the raw copy number as `seg.mean`, while `meta_seg.txt` and the docstring say `log2(CN/2)`.
- `data_cna.txt` is long format (`Hugo_Symbol, Sample_Id, Value`); `meta_cna.txt` uses `datatype: DISCRETE_LONG`.
- All deanon scripts warn on unmatched IDs but leave them unchanged.
- `vcf_to_seg.py` sets `num.mark=1` (ampliseq VCFs carry no probe-count).
- `meta_seg.txt`: `datatype: SEG`, `show_profile_in_analysis_tab: false`.

## Container Labels

Two process labels control container assignment in `nextflow.config`:

- `python` → `params.python_sif` (local Apptainer image built from `containers/python-ampliseq.def`)
- `mafsmith` → `params.mafsmith_container` (mafsmith + fastVEP, built from `containers/mafsmith-fastvep_v0.1.0-0.4.0.def`)

`params.mafsmith_data` (optional) is the mafsmith home with `GRCh37/{reference.fa,genes.gff3.gz}`. If it is null, empty or the string `"null"`, `DOWNLOAD_MAFSMITH` fetches it once into `assets/mafsmith` (storeDir). Either bundle then goes through **`PREPARE_MAFSMITH`** once per run, before any `MAFSMITH` task, which completes it in place: (1) it writes `GRCh37/reference.fa.fai` with an awk faidx, because `mafsmith fetch` never writes one and `mafsmith vcf2maf` fails without it (`Cannot read FASTA index`); (2) it decompresses `genes.gff3.gz` → `genes.gff3` atomically, and redoes it when the size doesn't match the `.gz`. mafsmith decompresses on first use only if `genes.gff3` is absent, non-atomically, so parallel tasks used to read a half-written GFF and silently lose annotation (`-`/`IGR`, no error). The rewrite's new mtime makes fastVEP rebuild `genes.gff3.fastvep.cache`, whose writes are atomic and checked by mtime. Both steps write into the bundle, so it must be writable. `MAFSMITH` fails early if `reference.fa`, `reference.fa.fai` or `genes.gff3.gz` is missing or empty, unless `ext.args` contains `--skip-annotation` (the module test relies on this). **Gotcha:** `-stub-run` with `skip_vcf2maf=false` stores an _empty_ `assets/mafsmith`, which storeDir then reuses; delete it after stub runs. **Contig naming:** the bundle is assumed to use `mafsmith fetch` (Ensembl) naming: `1, 2, X, MT`. mafsmith refuses a VCF whose contigs are not in the FASTA (`VCF chromosome 'chr1' was not found in the FASTA index`); the GFF's naming does not matter (verified on real Ensembl GRCh37 chr12). `MAFSMITH` therefore strips `chr` from the PASS VCF (records + `##contig`, `chrM`→`MT`) and puts it back on the MAF `Chromosome` column (`MT`→`chrM`) when the source VCF was chr-prefixed, as Pisces VCFs are. A pre-staged `--mafsmith_data` with chr-named FASTA is not supported. `Entrez_Gene_Id` is always 0 with mafsmith/fastVEP (no Entrez source). **Nirvana CSQ fields:** ampliseq/Pisces VCFs come Nirvana-annotated with INFO `CSQR`/`CSQT`. mafsmith 0.1.0 finds fastVEP's header with `contains("ID=CSQ")` (first match wins), so it parsed fastVEP's annotations with CSQR's 3-field layout and every row came out `Hugo_Symbol=Unknown` / `Targeted_Region`. `MAFSMITH` therefore drops all `CSQ*` INFO headers and values from the PASS VCF. fastVEP writes a fresh `CSQ` of its own, so nothing is lost; the `SAMPLE_NIRVANA` fixture covers this. Otherwise, all-`Unknown`/`Targeted_Region` output means `--skip-annotation` (no CSQ header). With annotation on, a record without CSQ is dropped, an intergenic one gets `-`/`IGR`, and a transcript without SYMBOL gets its `ENST` ID. Incremental runs never re-annotate: samples whose four per-sample files exist in `--outdir` are skipped and their old `_mutations.txt` merged as-is. Mutation_Status stays blank (tumor-only). `MERGE_MUTATIONS` merges per-sample files by column name, so older vcf2maf-era files (different columns) still merge correctly.

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

**Filter mode is applied per sample**, so switching `--filter_tsv_variants` or changing a threshold does not touch samples already in `--outdir`: use a fresh `--outdir` (or delete `samples/<id>/`) to re-filter them.

---

## Standalone Scripts

All Python scripts write output relative to `os.getcwd()` — run from the target output directory:

```bash
cd /path/to/output
python3 /path/to/bin/format_tsv_to_sv.py    <export.tsv>  <SAMPLE_ID>
python3 /path/to/bin/format_cna.py    <export.tsv>  <SAMPLE_ID>
python3 /path/to/bin/vcf_to_seg.py    <cnv.vcf>     <SAMPLE_ID>
python3 /path/to/bin/format_mutations.py data_mutations.txt <linking_file>
python3 /path/to/bin/format_sv.py       data_sv.txt         <linking_file>
python3 /path/to/bin/format_cna_deanon.py data_cna.txt      <linking_file>
python3 /path/to/bin/seg_deanon.py     data_seg.txt         <linking_file>
```

`bin/run_pipeline.sh` orchestrates all of the above — contains **hardcoded cluster paths** that must be updated before use.
