# CRCHUM-CITADEL/nextflow-cbioportal (oncoanalyser): Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Added

- **SAGE/PAVE INFO annotations are carried into the MAF.** `mafsmith --retain-ann` only reads CSQ subfields, so every plain INFO field the caller wrote was silently dropped — including PAVE's gnomAD frequency, which was sitting in the input VCF as `GND_FREQ` the whole time. A new `bin/annotate_maf_with_vcf_info.py` step in `MAFSMITH` joins ten of them back: `GND_FREQ` as **`gnomAD_AF`** (the conventional MAF spelling) plus `HMF.TIER`, `HMF.CLNSIG`, `HMF.CLNSIGCONF`, `HMF.PON_COUNT`, `HMF.MAPPABILITY`, `HMF.MSG`, `HMF.TNC`, `HMF.REP_C` and `HMF.MH`. `meta_sequenced.txt` now declares `namespaces: HMF` so cBioPortal surfaces the prefixed columns instead of ignoring them.
  - The join is keyed on `(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2)` with the VCF side normalised exactly the way mafsmith does it (`src/vcf/normalization.rs`: strip the common REF/ALT prefix advancing POS, render an emptied allele as `-`, then `Start = pos-1` for an insertion and `pos` otherwise). Keying on raw VCF `POS`/`REF`/`ALT` would have matched SNVs only and silently skipped every indel.
  - The columns are always written, even when a VCF declares none of the fields, so every per-subject MAF keeps identical width and the group-level `collectFile` merge cannot go ragged.
  - The step reports its match rate and **fails** when no MAF row matches any VCF record, since a silently 0% join means the key is wrong. `--min-match-rate` can raise that floor.
  - **Incremental runs need attention.** A subject cached from a previous run keeps its older, narrower MAF, and mixing those with freshly annotated ones in the group merge produces a ragged `data_mutations_dna_rna_germline.txt`. Delete the per-subject output directories and reprocess, exactly as the 4.0.0 chrX/chrY change required.
  - Covered by `tests/test_annotate_maf_with_vcf_info.py` (11 cases: SNV/insertion/deletion coordinates, multi-allelic, `chr`-prefix mismatch, width invariants, and both failure modes).

### Fixed

- **Somatic mutations showed ~0% allele frequency in cBioPortal.** cBioPortal derives allele frequency from `t_alt_count / (t_alt_count + t_ref_count)`, and since the vcf2maf → mafsmith switch those columns held the **normal** sample's counts. PAVE somatic VCFs list the normal column first and the tumor second. vcf2maf.pl defaulted `--vcf-tumor-id` to `--tumor-id`, so passing the barcodes also picked the columns; mafsmith keeps the two separate, and without `--vcf-tumor-id` it takes the first sample column as the tumor. `t_*` and `n_*` were therefore swapped on every row — `t_alt_count` was 0 for somatic calls. `MAFSMITH` now passes `--vcf-tumor-id` / `--vcf-normal-id` alongside the barcodes. The swap also fed `dna_vaf` into the ML vaf/integrated/hybrid mutation encodings and mafsmith's choice of the somatic ALT allele, which are corrected by the same change. Covered by the new `tests/modules/mafsmith.nf.test`, the first test to run the real mafsmith binary (with `--skip-annotation`, via a new `task.ext.args` hook), which fails against the previous behaviour. **Cached subjects keep their swapped MAFs** — `resolveSubjectCache()` matches on file existence, not content — so delete the per-subject output directories and reprocess.
- **`container_sigprofiler` pointed at a package that does not exist.** The default was `oras://ghcr.io/crchum-citadel/sdp-sigprofiler:1.1.3`, but the image is published as `crchum-citadel/sigprofiler:1.1.3` (`sdp-sigprofiler` returns 404). Any run without a pre-cached image failed to pull it. The default and the `test` profile now use the published name.
- **A gene with an unparseable `minCopyNumber` was written as a high-level amplification.** `gen_purple_cnv_to_cbioportal.R` coerces the column with `as.numeric()`, which yields `NA` on anything malformed, and then discretises with data.table's `fcase()`. `fcase()` returns its `default` when every condition is `NA` rather than propagating `NA`, so the gene fell through to `Value = 2` — byte-identical to a real high-level amplification. Such genes are now dropped, with a counted `warning()` naming how many. They are deliberately **not** written as `NA`: `NA` passes `validateData.py` (it is listed in `CNADiscreteLongValidator.ALLOWED_CNA_VALUES`) but aborts the load, because `CnaUtil.createAlteration()` ends in `Integer.valueOf(value).shortValue()` with no `NA` branch and the line loop in `ImportCnaDiscreteLongData` has no `try`/`catch`. An absent (gene, sample) pair is the format's own way of saying "not profiled" — the importer folds `DISCRETE_LONG` into the wide `DISCRETE` form and renders a missing pair as an empty cell. Covered by `tests/modules/purple_cnv_to_cbioportal.nf.test`, which fails against the previous behaviour. Output is unchanged for any cohort whose PURPLE files parse cleanly, so no checksum moves in the test fixtures

## 4.0.0 - 2026-09-17

### Breaking

- **Supports oncoanalyser 3.0 only.** Output from oncoanalyser 2.3 is no longer accepted: the Isofox fusion converter hard-requires the 3.0 `pass_fusions.tsv` columns (`Name`, `TranscriptUp/Down`, `ExonUp/Down`, `SplitFrags`/`RealignedFrags`/`DiscordantFrags`) and stops on anything else, and the genomic workflow resolves the 3.0 directory layout (`pave/`, `sage_append/<sample>-T/`, `.esvee.somatic.vcf.gz`, `.isf.*.tsv`). No parameter gates this — re-run oncoanalyser 3.0 rather than downgrading the pipeline.
- **Requires Nextflow >= 26.04.4** (manifest `nextflowVersion`, CI, and the documented `module load`).
- **Site-specific defaults have moved out of `nextflow.config`** into `nextflow_citadel.config`, loaded by the new `citadel` profile. CRCHUM runs must add it — `-profile slurm,apptainer,citadel`, listed last. `nextflow.config` now ships portable defaults only: public `oras://` container images, empty reference-data paths, `outdir = 'output'`, no notification address. `ensembl_annotations`, `genome_reference` and `mafsmith_data` therefore have no default and must be supplied for genomic/both mode.
- **`vcf2maf` is replaced by `mafsmith`** ([nf-osi/mafsmith](https://github.com/nf-osi/mafsmith) v0.1.0). The `VCF2MAF` local module, the nf-core `vcf2maf` module and the vcf2maf container definition are gone, along with the params `container_vcf2maf`, `vep_params` and `vep_path`, which nothing else read. Two new params replace them: `container_mafsmith` and `mafsmith_data`, the mafsmith home directory holding its reference data (gene GFF3 + reference FASTA), symlinked in as `.mafsmith`. `mafsmith_data` is optional and auto-downloads like `vep_data` and `pcgr_data` (see Added), though pre-staging is preferred because its reference must match the genome SAGE called against. `vep_data` stays: PCGR still uses it, and `genome_reference` is still required, now only by SigProfiler. **`data_mutations_dna_rna_germline.txt` gets narrower:** mafsmith emits 53 MAF columns against vcf2maf's ~133, so VEP passthrough fields such as `Consequence`, `IMPACT`, `SYMBOL`, `CCDS` and `RefSeq` are no longer written. Nothing in the pipeline reads them — both MAF-consuming R scripts address columns by name and use none of the missing ones — but cBioPortal displays some. The new `mafsmith_retain_ann` param puts the useful ones back via `mafsmith vcf2maf --retain-ann`, defaulting to `Consequence,IMPACT,SYMBOL,Feature_type,CCDS,EXON,INTRON,Existing_variation,CLIN_SIG,AF` — every name taken from fastVEP's `DEFAULT_CSQ_FIELDS`, so none of them is a guess. Set it to `""` for the bare 53. vcf2maf's `RefSeq` and `VARIANT_CLASS` have no fastVEP counterpart and cannot be restored this way.
- **`ensembl_annotations_expr` is removed.** Isofox expression now maps Ensembl→Entrez through `ensembl_annotations`, the same BioMart TSV the CNV and SV converters already used, so only one annotation file has to be staged.
- **`generate_cancer_type` is replaced by `incremental`, with the opposite sense.** `cancer_type.txt` and `meta_cancer_type.txt` are now written on every run and _skipped_ by `--incremental`, which marks a follow-up load into a study cBioPortal already holds (re-registering the cancer type errors there). A first load needs no flag; per-subject output caching remains automatic and unrelated.
- **CNV copy-number calls on chrX/chrY are now scored against the sample's own sex, not a flat diploid baseline.** PURPLE reports absolute copy number, which is haploid on chrX/chrY in a male sample — the converter previously divided everything by 2 regardless, so every male sample's chrX and chrY came out as a whole-chromosome hemizygous loss (`seg.mean -1.0`, discrete call `-1`), and every female sample's chrY came out as a homozygous deletion of every gene on it (`seg.mean -10.9658`, discrete call `-2`). `gen_purple_cnv_to_cbioportal.R` now reads `gender` from the sample's `purple/<subject>-T.purple.purity.tsv` (PURPLE's own MALE/FEMALE/MALE_KLINEFELTER call) and scores chrX/chrY against the correct expected copy number (1 for a male's X/Y, 2 otherwise); chrY rows are dropped entirely for a female sample rather than scored as a deletion. Autosome output is unchanged. A missing or unreadable purity file falls back to the old diploid-everywhere behaviour with a warning, so this is a no-op for a subject with no purity file. **`data_cna_hg38.seg` and `data_cna_long.txt` checksums change for every cohort containing a male sample or a female sample with a genotyped chrY.** Incremental processing does not pick this up automatically — `resolveSubjectCache()` matches on file existence, not content, so already-published subjects must have their per-subject output directory deleted and be reprocessed to get corrected chrX/chrY calls.
- **Germline variants now actually reach `data_mutations_dna_rna_germline.txt`.** `CONVERT_CPSR_TO_MAF` ran `gen_convert_cpsr_to_maf.R` to write `tmp.<sample>.somatic_rna_germline.maf` and then immediately `mv`-ed the somatic+RNA input over that same path, discarding the merged result; every published `.somatic_rna_germline.maf` was therefore somatic+RNA only, despite its name. The `mv` is gone (`dragen_4.4`'s otherwise-identical module never had it). Expect the mutation file to grow by the CPSR Pathogenic / Likely_Pathogenic / VUS calls per sample, so its checksum changes.
- `data_clinical_sample.txt` gains two columns (`CANCER_TYPE_LABEL`, `TUMOR_TISSUE_SITE_LABEL`), so its checksum changes.

### Added

- `CANCER_TYPE_LABEL` and `TUMOR_TISSUE_SITE_LABEL` columns in `data_clinical_sample.txt`: human-readable labels for the ICD-O topography codes already carried by `CANCER_TYPE_CODE` (`cancer_type_code`) and `TUMOR_TISSUE_SITE` (`specimen_anatomic_location`), reverse-mapped through the MOHCCN primary-site table the same way `RELAPSE_SITE_LABEL` already was. Labels resolve to the C+2-digit parent code, since the MOHCCN table only carries parent codes (`C22.0` and `C22.1` both label as "Liver and intrahepatic bile ducts")
- Timeline generation from ARGO clinical CSVs: produces a single combined `data_timeline.txt` with all event types (surgery, treatment, status, specimen, lab*test), built by the `GENERATE_TIMELINE*\*`modules and unioned by`MERGE_TIMELINE`
- MOHCCN mapping tables (`assets/mohccn_*` CSVs) for translating plain-text values to ICD-O / ontology codes in clinical and timeline output
- New params: `mohccn_primary_site_map`, `mohccn_specimen_tissue_source_map`, `mohccn_treatment_intent_map`
- MOHCCN code integration in the clinical output (primary site, tissue source, treatment intent)
- Specimen SAMPLE_TYPE derivation from nearest follow-up disease status (Primary / Recurrence / Metastasis)
- Biomarker lab test pivoting to cBioPortal long format (PSA, CEA, CA125, ER/PR/HER2, HPV)
- DNA vs RNA SV distinction (`ecb509b`): DNA_Support/RNA_Support flags properly set per source
- `combine_cbioportal_outputs.py`: support for `data_timeline.txt` (union-of-columns append), `cancer_type.txt`, per-subject folders, and glob-discovered `meta_*.txt` / `case_lists/*.txt` so new outputs are never silently dropped
- `combine_cbioportal_outputs.py`: `--strict` flag, and a warning listing any input file no merge strategy handles
- `tests/test_combine_cbioportal_outputs.py` pytest suite, plus a pytest job in the linting workflow
- Per-module nf-tests for every clinical module (there were none): `build_clinical_table`, `write_clinical_{sample,patient}`, `generate_timeline_{surgery,treatment,status,specimen,lab_test}` and `merge_timeline`, including empty-input and deterministic-ordering cases
- `citadel` profile and `nextflow_citadel.config` for CRCHUM site settings; `SITE CONFIGURATION` section in `docs/usage.txt` for everyone else
- `incremental`, `mafsmith_data` and `container_mafsmith` declared in `nextflow_schema.json` (an undeclared param trips nf-schema validation)
- `DOWNLOAD_MAFSMITH`: when `mafsmith_data` is empty the pipeline provisions it with `mafsmith fetch --genome grch38 --ensembl-release 113`, stored in `assets/mafsmith/` and gated on there being a subject to process, exactly like the VEP and PCGR downloads. **Pre-staging is still preferred:** the fetch pulls Ensembl's primary assembly, whose contigs are `1`/`2`/`X`, while oncoanalyser aligns against GATK hg38 (`chr1`/`chr2`/`chrX`), so the fallback suits testing and portability rather than a real cohort. `docs/usage.txt` documents building a bundle by hand
- `fastvep` is built into the mafsmith container. `mafsmith vcf2maf` falls back to `which fastvep` when `<data-dir>/bin/fastvep` is absent, so `DOWNLOAD_MAFSMITH` passes `--skip-fastvep` and no cargo toolchain is needed at run time. **`containers/mafsmith_v0.1.0.def` must be rebuilt and repushed for this**

### Changed

- Clinical output is generated by focused modules instead of two monolithic ones, so a failure identifies the output it was building. `FORMAT_CLINICAL` (which ran twice per group, repeating ~95% of its work per mode) becomes `BUILD_CLINICAL_TABLE` + `WRITE_CLINICAL_SAMPLE`/`WRITE_CLINICAL_PATIENT`; `GENERATE_TIMELINE` becomes five per-EVENT_TYPE modules plus `MERGE_TIMELINE`. Output is unchanged.
- Shared R helpers extracted to `bin/clinical_common.R`, staged into each process as an explicit `path()` input so it participates in the task hash
- Removed dead `modules/local/assign_date` and `modules/local/concat_results`
- Renamed the stale `CRCHUM-CITADEL/nextflow-sante-precision` to `nextflow-cbioportal` throughout
- The pipeline directory is now named for the upstream tool version it consumes: `oncoanalyser/` → `oncoanalyser_3.0.0/` (its sibling `dragen/` → `dragen_4.4/`). CI working directories and the monorepo README/CLAUDE.md tables follow
- `process_medium_memory` drops from 60 GB / 36 h to 36 GB / 23 h, which is what the mutation step actually needs now that mafsmith has replaced vcf2maf + VEP
- `tests/subworkflows/genomic_mutations.nf.test` runs under `-stub-run` and asserts the wiring instead of MAF content, with a second case covering the download fallbacks
- `MAFSMITH` locates the VCF sample columns by the `FORMAT` header and the MAF `Tumor_Sample_Barcode` column by name, instead of hardcoding VCF fields 10/11 and MAF field 16. Verified to select the same columns on the existing fixtures, and it now fails loudly if the MAF header lacks the column
- `gen_convert_cpsr_to_maf.R` wrote germline rows wider than the header they were written under. It builds each row as a list keyed on the incoming MAF's columns, but `maf_entry$X <- v` appends when `X` is not already a slot, so the six assignments naming `Feature_type`, `Consequence`, `IMPACT`, `FILTER`, `CCDS` and `RefSeq` widened the row whenever the MAF did not already carry them — and `rbind()` onto a zero-row frame adopts the wider shape rather than erroring. vcf2maf's ~133 columns happened to cover all six; mafsmith emits 53, which turned this into ragged output (59 fields against a 53-field header, reproduced in the R container). Each row is now pinned back to the header's columns before it is appended

### Fixed

- `optional: true` on the timeline outputs was written as an argument to `path()` rather than an option on the output declaration, so it had no effect: a process whose event type produced no rows failed with `MissingFileException` despite exiting 0. The pre-existing `GENERATE_TIMELINE` had the same latent mistake
- `data_timeline.txt` row order is now deterministic. `MERGE_TIMELINE` stable-sorts by a fixed EVENT_TYPE order, since `mix()`/`groupTuple()` guarantees no arrival order across the five per-event channels
- `MAFSMITH` has a `stub:` block. A process without one runs its **real** script under `-stub-run` instead of erroring, which would have broken `tests/incremental_pipeline.nf.test`
- The empty `params { }` block in `tests/nextflow_subworkflow.config` was dead, and Nextflow rejects it (`Unknown config attribute 'params'`) as soon as a test adds a config of its own
- `workflows/genomic.nf`'s path-convention comment described the pre-3.0 layout while the code below it read the 3.0 one
- `<study>.tar.gz` in `both` mode was genomic-only: `PACKAGE_CBIOPORTAL` now runs from `main.nf` after `CLINICAL`, so the archive also carries the clinical files, `data_timeline.txt`, `meta_timeline.txt` and `util_linking_file.txt`
- `data_mutational_signatures_counts_ID.txt` was published to the output directory but missing from `<study>.tar.gz`
- `combine_cbioportal_outputs.py`: clinical and timeline merges no longer misalign rows when the two batches carry different column sets
- `combine_cbioportal_outputs.py`: a genuine expression value of `0` is no longer replaced by the fill value in wide-matrix merges
- `combine_cbioportal_outputs.py`: row-append no longer glues two rows together when the first file lacks a trailing newline (`util_linking_file.txt`)
- `combine_cbioportal_outputs.py`: in-place merge keeps a backup until the swap succeeds instead of deleting the target first
- `tests/.nftignore` excludes `*.tar.gz`: tar embeds per-file mtimes, so the archive md5 differed on every run
- Timeline START_DATE defaults to 0 (diagnosis day) when source date is missing
- Clinical sample deduplication: better merging to remove duplicate sample-level rows
- Removed `analyte_type` from sample-level clinical output rows
- Status data frame handling in clinical output
- MOHCCN map lookups now trim whitespace on both the label and code columns consistently between the clinical and timeline outputs (the clinical map previously did not trim, so a whitespace-padded map entry could resolve in `data_timeline.txt` but silently become `NA` in `data_clinical_sample.txt`); no checked-in `assets/mohccn_*` file is actually affected today, but the two lookups no longer risk drifting apart

## v1.0.0dev - [date]

Initial release, created with the [nf-core](https://nf-co.re/) template.

### Added

- Genomic mode: mutations (vcf2maf), CNV (PURPLE), SV (ESVEE + Isofox fusions), expression (Isofox), mutational signatures (SBS/DBS/ID)
- Clinical mode: ARGO/ICGC-ARGO CSV processing into cBioPortal patient/sample attribute files
- Both mode: automatic linking of genomic sample IDs to clinical patient IDs
- Incremental processing: skip already-processed subjects on re-runs
- ML feature table generation (CNV, expression, mutation, SV)
- Samplesheet generators for genomic and clinical input
- cBioPortal output packaging (.tar.gz)
