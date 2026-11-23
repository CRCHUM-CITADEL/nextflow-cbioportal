# CLAUDE.md — oncoanalyser_3.0.0

Oncoanalyser 3.0 WGS/WTS + clinical CSVs → cBioPortal. Nextflow >= 26.04.4. 2.3 output is
rejected (converters `stop()` without the 3.0 columns).

## Mutations (MAF)

- Inputs are PURPLE's VCFs (`purple/<subject>-T.purple.{somatic,germline}.vcf.gz`; PAVE-annotated, plus
  PURPLE INFO), also for SigProfiler DBS/ID. `pave/` is no longer read. RNA stays `sage_append/`
- VCF→MAF is `MAFSMITH` (mafsmith + fastVEP), not vcf2maf. mafsmith reads only `$HOME/.mafsmith`,
  so the module sets `HOME=./` and symlinks `mafsmith_data` there. `genome_reference` is now used
  only by SigProfiler
- `mafsmith_data` comes from the resource kit. `DOWNLOAD_MAFSMITH` (`mafsmith fetch`) is test-only:
  it gets Ensembl contigs (`1`), not GATK hg38 `chr1`. A hand-built bundle needs `--gff3` + `--ref-fasta`
- Production uses the kit's `bin/fastvep` (`<data-dir>/bin/fastvep` wins over the container's 0.4.0): the
  one `mafsmith fetch` installs (`cargo install` from fastVEP `master`; 0.5.0 in kit 1.1.0; mafsmith 0.1.0
  is still the latest). Changing it changes output and needs a new kit version. fastVEP writes
  `genes.gff3` + `genes.gff3.fastvep.cache` beside the GFF3 on first use, so the kit ships them pre-built.
  Image def: `containers/mafsmith-fastvep_v0.1.0-0.4.0.def`, tag `mafsmith-fastvep:<mafsmith>-<fastvep>`
- mafsmith does not derive `--vcf-tumor-id` from `--tumor-id`: without `--vcf-tumor-id` /
  `--vcf-normal-id` it takes the first VCF sample (the normal) as tumour and `t_*`/`n_*` swap.
  `MAFSMITH` passes both; it finds sample columns via the `FORMAT` header and MAF columns by name.
  `mafsmith.nf.test` runs the real binary with `--skip-annotation` to guard this
- `Mutation_Status` is only `Somatic` (set by `MAFSMITH`) or `Germline` (CPSR rows); cBioPortal and
  `ml_format_mutation.R` rely on it. RNA evidence goes in `RNA.depth/vaf/ref_count/alt_count/expressed`
  (namespace `RNA` in `meta_sequenced.txt`), never in `Mutation_Status`. Plus a summary column
  `RNA_expressed` = `RNA(d=<depth>,v=<vaf>)` / `RNA(data=NO)` (empty on germline rows)
- Every DNA subject reaches the MAF with or without RNA (`remainder: true` join; `[]` RNA input);
  the `RNA.*` columns are always written so the group merge stays rectangular
- `--retain-ann` reads CSQ subfields only, so `annotate_maf_with_vcf_info.py` joins 16 INFO fields
  back: `GND_FREQ` → `gnomAD_AF`, plus `HMF.{TIER,CLNSIG,CLNSIGCONF,PON_COUNT,MAPPABILITY,MSG,TNC,REP_C,MH}`
  and PURPLE's `HMF.{PURPLE_AF,PURPLE_CN,PURPLE_VCN,PURPLE_MACN,BIALLELIC,SUBCL}` (a flag reads `1`).
  Its join key must normalise alleles exactly like mafsmith (`src/vcf/normalization.rs`: strip common
  prefix advancing POS, empty allele → `-`, insertion `Start = pos-1`) or every indel is silently missed.
  Columns are written unconditionally; the step fails when nothing matches
- mafsmith writes `Entrez_Gene_Id = 0`; `fill_maf_entrez.py` fills it from `ensembl_annotations`
  (Ensembl `Gene` id, then `Hugo_Symbol`; ambiguous keys stay `0`). CPSR rows use CPSR's `ENTREZ_ID`
- mafsmith writes 53 columns (vcf2maf ~133). `mafsmith_retain_ann` names must match fastVEP's
  `DEFAULT_CSQ_FIELDS` (unknown → silently empty column); `RefSeq` / `VARIANT_CLASS` can't be restored

## Germline (PURPLE → PCGR/CPSR → CONVERT_CPSR_TO_MAF)

- PURPLE's germline VCF goes to CPSR **unfiltered**: no FILTER selection, no GT inference (the old
  `FILTER_GERMLINE_DNA` step was for PAVE, which wrote every GT as `./.`). CPSR itself keeps PASS only and
  drops `./.`, so PURPLE's `LOW_TUMOR_VCN` calls and its uncalled GTs never reach the MAF
- `PCGR` keeps only the normal's column, found by name (`<subject>-N`; PURPLE's column order is undocumented and
  can put the tumour first), and strips all INFO. No `<subject>-N` column fails the task
- A subject without a germline VCF also loses its somatic MAF (`MAFSMITH` input joins on `ger_dna_vcf`)
- INFO is stripped because CPSR aborts on clashing tags such as PAVE's `IMPACT`, kept by PURPLE. CPSR exits 0
  on failure, so any `- ERROR -` log line fails the task — never fall through to the header-only placeholder
- `PCGR` writes `<sample>.cpsr_debug/` (log, GT summaries, side outputs) in its work dir, **not published**: its copy
  hung shutdown on the cluster FS; `ext.args = '--debug'` keeps
  intermediates. CPSR output is named `<subject>-N`; `PCGR` renames it to `<sample>`, and writes a
  header-only placeholder only when CPSR wrote nothing (= no germline calls downstream)
- `CONVERT_CPSR_TO_MAF` keeps Pathogenic / Likely Pathogenic / VUS on CPSR's **final** classification
  (`CLASSIFICATION` in cpsr >= 2.3, `FINAL_CLASSIFICATION` <= 2.2, found by name, fails if neither),
  never `CLINVAR_CLASSIFICATION` (empty for novel variants). Accepts `Likely_Pathogenic` and
  `Likely Pathogenic`. Then drops `Intron`/`IGR`
- `gen_convert_cpsr_to_maf.R` pins germline rows to the incoming MAF's header: `maf_entry$X <- v` adds a
  slot for an absent `X` and `rbind()` widens silently. Add a column only if the annotator emits it
- PCGR 2.4.x needs bundle 20261106 + VEP 116 — bump `container_pcgr`, bundle and VEP cache together.
  Somatic (mafsmith GFF3), germline (VEP) and `ensembl_annotations` are all Ensembl 116: move them
  together (kit table: `bin/build_ensembl_annotations.py`, recipe in docs/resources.md)
- `PCGR` passes `--no_pdf` and `--clinvar_cancer_only` (2.4-only flags), so every `container_pcgr`,
  `test` included, is 2.4.0. `fake_cpsr` logs ERROR without them; nf-tests never run real CPSR
- `DOWNLOAD_*` `storeDir`s carry the release (`assets/pcgr/20261106`, `assets/mafsmith/ensembl_116`) so an
  older download is never reused. `DOWNLOAD_VEP_TEST` is still the VEP 112 chr21 cache: wiring only

## SV / CNA / ML

- `data_sv.txt` needs Hugo symbols at both sites — unannotated rows are dropped
- DNA SVs come from `purple/<subject>-T.purple.sv.vcf.gz` (`gen_purple_sv_to_cbioportal.R`), not `esvee/`. PURPLE
  keeps every input record and adds `INFERRED` breakends and `TINC_RECOVERED` SVs; we keep `PASS` BNDs with a mate.
  The tumour column is the one named `<subject>-T` (column order is undocumented). Counts: `VF` → `Tumor_Variant_Count`,
  `SF` → split, `DF` → paired-end (ESVEE writes no `BVF`/`SR`/`PR`)
- BND strand: `(+,-)` DEL, `(-,+)` DUP, `(+,+)/(-,-)` INV, other chr TRANSLOC. The ALT-bracket regex needs
  `perl = TRUE` (in POSIX `\\[` inside `[]` is a literal backslash: every SV was INV). `DNA_Support=Yes, RNA_Support=No`.
  Isofox fusions: `Class=FUSION, DNA_Support=No, RNA_Support=Yes`
- `PURPLE_AF/CN/CN_CHANGE` are `[this breakend, mate]` on every record (PURPLE reverses the mate's), so each site
  takes element 1 of its own record → `HMF.PURPLE_{AF,CN,CN_CHANGE}_SITE{1,2}`, plus `HMF.PURPLE_JCN`. Every SV
  writer (Isofox too, as `NA`) emits all 7, or `MERGE_SAMPLE_SV` / the group `collectFile` go ragged.
  `meta_sv.txt` declares `namespaces: HMF`
- Isofox 3.0 `pass_fusions.tsv`: `Name` = `<up>_<down>` (split on first `_`); `Tumor_Variant_Count` =
  `SplitFrags + RealignedFrags + DiscordantFrags`. `Site1/2_Region`: transcript + exon → `Exon`, gene
  without exon ref → `Intron`, no gene → `NA` (UTR/Promoter not derivable)
- `gen_purple_cnv_to_cbioportal.R` scores against the expected baseline: autosomes 2, chrX 1 in a male,
  chrY 1 (2 if sex unresolved), chrY dropped in a female. Sex = `gender` in optional
  `purple/<subject>-T.purple.purity.tsv`; missing → diploid everywhere
- A gene whose `minCopyNumber` doesn't parse is **dropped**: `fcase()` defaulted it to `Value = 2`, and
  `NA` passes `validateData.py` but aborts the cBioPortal load. An absent (gene, sample) = not profiled
- `ml_format_cnv.R` / `ml_format_expression.R` require inputs named `data_cna_long.txt` / `data_expression.txt`
- `ml_mutation_processor.R` encodings are vectorized (the old row loop was quadratic). The vaf/effect max
  must propagate `NA` (no `na.rm`), like `max(current, NA)`. `fetch_hotspots()` uses
  `cancer_hotspots_data` if staged (`[]` when unset → live GET)
- `GENOMIC_ML` is opt-in (`run_ml`, default false). `both.nf.test` runs without it and asserts no `_ML_`
  task; ML is covered by `tests/subworkflows/genomic_ml.nf.test` only
- `cosmic_data` and `chimer_data` are all-or-nothing (empty `path` → `Path must not be empty`).
  ChimerKB is in the kit, so with `--run_ml` setting `cosmic_data` enables ML SV.
  `PIPELINE_INITIALISATION` rejects half-configuration; only `genomic_ml.nf.test` covers the branch

## Configuration

- Containers: R → `container_r`, Python → `container_python`, SigProfiler → `container_sigprofiler`,
  VCF→MAF → `container_mafsmith`. Compute nodes are offline (`NXF_OFFLINE=true`; pre-pull images)
- `nextflow.config` = portable defaults; CRCHUM values in `nextflow_citadel.config`
  (`-profile slurm,apptainer,citadel`, citadel last). `.gitignore` ignores `/nextflow_*.config` except it
- Reference data = one versioned kit (`--resources_dir`, docs/resources.md). Kit-backed params are not
  declared in `nextflow.config` (hidden in schema). Read them via `resource(name)`
  (`subworkflows/local/utils/main.nf`: explicit param, else kit file, else `""`), never `params.X`;
  `resourceLayout()` is the path table. Resolved at run time so profile values count. Run stops if
  `resources_dir/VERSION` != `params.resources_version`. `cosmic_data` is never in the kit. `test`
  uses the mini-kit in `assets/test_data/resources/`
- Optional params with a schema `format`/`pattern`/`exists` default to `null`, never `""` (nf-schema
  > = 2.7 rejects empty strings); test configs use `null` for "unset"; read by truthiness
- `--help` is nf-schema config help; `main.nf` returns on `helpRequested()` before initialisation
- `--incremental` = load into an existing study: skips `meta_study.txt`, `cancer_type.txt`,
  `meta_cancer_type.txt`
- New `bin/` scripts need `chmod +x` (else exit 126). Call them by bare name; the shebang picks the interpreter

## Clinical

- `bin/clinical_common.R` holds shared helpers and is staged as an explicit `path()` input so it is in the task hash
- `BUILD_CLINICAL_TABLE` (once per group, writes `clinical_merged.tsv` with `na="NA"`, read back as
  character) → `WRITE_CLINICAL_SAMPLE` / `WRITE_CLINICAL_PATIENT` apply `col_defs` only
- `GENERATE_TIMELINE_{SURGERY,TREATMENT,STATUS,SPECIMEN,LAB_TEST}` → `MERGE_TIMELINE` → `data_timeline.txt`.
  Parts write nothing when empty. `MERGE_TIMELINE` owns column order (PATIENT_ID, START_DATE, STOP_DATE,
  EVENT_TYPE, then alphabetical — contract with `combine_cbioportal_outputs.py`) and sorts rows by that
  EVENT_TYPE order for determinism. Missing START_DATE → 0
- `is_cbio_sample()` is the only sample-eligibility rule (`Total DNA` + `Tumour` + `Solid tissue`,
  `-\d+[A-Z]*[DR]T$`): normals never become sample rows; donors without a tumour drop out. Timeline
  events are filtered by patient only; `genomic_subjects` restricts to genomic subjects. Nearest
  follow-up sets specimen SAMPLE_TYPE
- Each sample row takes the diagnosis of its own specimen (`submitter_primary_diagnosis_id`), else the
  patient's first. Other tables merge by patient alone
- MOHCCN maps (`mohccn_*_{primary_site,specimen_tissue_source,treatment_intent}.csv`) go through
  `read_mohccn_map()` / `reverse_mohccn_map()`. `apply_ps_label()` adds `*_LABEL` to `RELAPSE_SITE`,
  `CANCER_TYPE_CODE`, `TUMOR_TISSUE_SITE`, at C+2-digit parent level, joining multi-site lists with `"; "`.
  `_LABEL` avoids cBioPortal's reserved OncoTree `CANCER_TYPE*` columns
- Dashboard priority (last `col_defs` field): `PRIMARY_SITE` 1500 (must stay > 1000, the Data Types
  chart), others 1; `0` hides. `write_cbio_table()` silently drops a col_def whose source is missing —
  `write_clinical_sample.nf.test` pins priorities and labels
- `--mode clinical` needs `--study_id` and `--clinical_samplesheet` or `--sample_registrations`; it makes
  no archive / `meta_study.txt` / case lists. An optional `--genomic_samplesheet` filters to its subjects
  (TSV shaped like `util_linking_file.txt`, kept in work dir). Clinical test configs must set
  `genomic_samplesheet = null` (the `test` profile sets one)

## Re-running (`-resume`)

- Nextflow's `-resume` is the only caching: every samplesheet subject feeds the tasks on every run, and
  unchanged tasks come back from `work/`. There is no skip logic reading `outdir`
- Task inputs must be deterministic or `-resume` reruns them: after a `groupTuple()` pass a meta built from
  stable fields (`[group: group]` for group merges), never `metas[0]` (arrival order), and sort file lists.
  Only `PACKAGE_CBIOPORTAL` reruns (published with `mode: 'move'`); `tests/resume.nf.test` guards this
- The samplesheet defines the cohort: group files are rebuilt each run; a dropped subject is warned about
  (via the previous `util_linking_file.txt`) and left out of the study
- No somatic VCF → no reference download (`DOWNLOAD_*` gated on `count().filter { it > 0 }`)
- Case lists (`case_list_defs` in `GENOMIC_AGGREGATE_OUTPUT`): standard suffixes `_sequenced`, `_cna`, `_cnaseq`, `_sv`,
  `_rna_seq_v2_mrna`, `_3way_complete`, each with its `case_list_category`; combined lists are intersections. Built from one
  `[profile, sample]` stream via `toList()` (not `collect()`: must emit when empty). Empty lists are never written. No
  `cases_all.txt`: `add_global_case_list` makes `_all`. On `--incremental`, cBioPortal drops loaded samples from every list
  not in `case_lists/`, so a list must hold exactly the samples with that data. RNA = samples in `data_expression.txt`

## `--mode genomic-dna-only`

- Modes are defined only in `pipelineModes()` / `runsGenomic()` / `isDnaOnly()` (`subworkflows/local/utils/main.nf`);
  never compare `params.mode` to a genomic list directly. Matched exactly: CLI params are immutable
- RNA channels (`sage_append`, Isofox) are `channel.empty()`; the MAF still gets empty `RNA.*` columns.
  `GENOMIC_AGGREGATE_OUTPUT` drops TPM input and the expression meta file
- `hasOncoOutput(meta, rel)` (`workflows/genomic.nf`) skips a step when `purple/`, the PURPLE SV VCF or
  `sigs/` is absent
- Under dna-only only, meta files exist only for data types with output (labels travel with their
  lists; case lists do this in every mode), and groups come from CNV/SV/mutation output
- Every dna-only branch is behind `isDnaOnly()`: other modes must stay byte-identical
- `data_clinical_sample.txt` = TMB/HRD attributes, as in genomic mode. No DNA-only `both`

## Outputs

- `PACKAGE_CBIOPORTAL` is called from `main.nf` with `GENOMIC` + `CLINICAL` `package_files`, so in
  `both` mode the archive mirrors the loadable study. Its md5 isn't reproducible (`*.tar.gz` in
  `tests/.nftignore`) — assert on contents. Published with `mode: 'move'` (a large `copy` stalled
  shutdown on "Waiting files transfer"), so it re-runs on `-resume` and tests read it from `outdir`
- `bin/combine_cbioportal_outputs.py` (standalone) merges two study folders/archives per file category
  (row-append, MAF 2-line header, wide matrix, timeline/clinical/SV column union, `meta_*`/case lists by
  glob). Not `machine_learning/` (cohort-normalised). Unhandled files are reported (`--strict` fails);
  overlapping sample IDs are an error. New outputs → its dispatch table + `tests/test_combine_cbioportal_outputs.py`
- Clinical/genomic column labels and descriptions never name the tool (no "(PURPLE)", "(CHORD)")
- TMB/HRD sample attributes: `TMB_PER_MB, TMB_STATUS, TML, TML_STATUS` from `purple.purity.tsv`,
  `HRD_SCORE` (`hrd`, else `p_hrd`), `HR_STATUS, HRD_TYPE` from `chord/<subject>-T.chord.prediction.tsv`.
  Columns by name; missing file → `NA`, missing column → fail. `bin/gen_genomic_sample_attributes.py` owns
  the col defs; `GENOMIC_AGGREGATE_OUTPUT` merges them (`collectFile(skip: 5)`) into an intermediate, not a
  study file. **cBioPortal rejects a second `SAMPLE_ATTRIBUTES` file** ("Multiple sample attribute files"),
  so: genomic and genomic-dna-only modes publish it as `data_clinical_sample.txt` + meta (`workflows/genomic.nf`); both mode
  passes it to `CLINICAL`, whose `WRITE_CLINICAL_SAMPLE` / `GENERATE_CLINICAL_TEMPLATE` run
  `join_genomic_sample_attributes.R` (by `SAMPLE_ID`, before `SPLIT_CLINICAL`). Clinical mode passes `[]`
- Signatures: `docs/mutational_signatures.md`. Process labels: `conf/base.config`

## Testing

- `nf-test test tests/<x>.nf.test --profile=+apptainer`. Channel file outputs are `String`: use `path(f.toString())`
- Sort snapshots by basename; `collectFile` with `storeDir` won't create dirs — `mkdirs()` in setup
- A process without `stub:` runs its real script under `-stub-run` — every genomic-path module needs one
- Pipeline tests run every task for real except CPSR (`fake_cpsr` stand-in in the bcftools image, via
  `containerOptions --env PREPEND_PATH`: Nextflow starts containers with a clean env) and fastVEP
  (`mafsmith --skip-annotation`); see `tests/nextflow_both.config`. Wiring-only tests use
  `tests/nextflow_stub.config` under `-stub-run` (container engines off)
- `tests/resume.nf.test` runs the stub pipeline twice (setup + `-resume`, opposite task orders via
  `tests/nextflow_resume.config`) and fails on any rerun task but `PACKAGE_CBIOPORTAL`
- `WorkflowTask` has only `name` and `success`: `workflow.trace.succeeded().collect { it.name }`
- `${projectDir}` in a class-level `def` string inside a test isn't interpolated — inline it in the process block
- `bin/create_test_data.sh` must reproduce `assets/test_data/` byte for byte (no `git status` change after a
  run); edit a generated fixture there, not by hand. Layout: `assets/test_data/README.md`
- Test samplesheets with relative `folder`s don't resolve under nf-test — write an absolute one outside
  `outputDir` (it would enter the snapshot)
