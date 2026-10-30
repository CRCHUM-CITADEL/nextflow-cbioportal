# crchum-citadel/ampliseq-cbioportal

## Introduction

**crchum-citadel/ampliseq-cbioportal** formats ampliseq genomic data (VCFs + structural variant TSVs) into cBioPortal-compatible files: mutations (MAF), copy number alterations, structural variants, segmentation data, and clinical data.

Requires Nextflow >= 26.04.4 and Apptainer.

## Usage

### 1. Prepare input files

#### Generating a samplesheet

A helper script can automatically generate the samplesheet from a directory of sample folders:

```bash
python3 bin/generate_samplesheet.py /path/to/sample_data -o samplesheet.csv
```

The input directory can contain either cohort subdirectories (each with per-sample folders) or sample folders directly. Each sample folder must contain `*-basespace-pisces.final.vcf.gz` and `analysis_*_export.tsv`. The cohort subdirectory name becomes the `group` column, the `sample_id` is extracted from the VCF filename prefix (before `-basespace-`), and the `subject_id` is derived from the sample*id (part before the first `*`). Samples missing required files are skipped with warnings.

#### Samplesheet format

**Samplesheet** (`samplesheet.csv`):

```csv
group,subject_id,sample_id,folder_location
cohort_A,PATIENT_001,SAMPLE_001,path/to/samples/SAMPLE_001
```

Each sample folder must contain:

- `analysis_*_export.tsv` — structural variant / CNA export
- `*-basespace-pisces.final.vcf.gz` — compressed VCF for mutation calling
- `*-basespace-cnv.final.vcf` — CNV VCF for segmentation (`CN` FORMAT field required)

**Linking file** (`linking_file.txt`, tab-separated) — maps anonymized → real IDs:

```
sample_id	deanon_sample_id	deanon_patient_id
SAMPLE_001	PATIENT_001	PATIENT_001
```

One patient may have multiple rows (one per sample). `deanon_patient_id` is used to filter `data_clinical_patient.txt` to only patients whose samples are in the samplesheet.

**Patient file** (tab-separated): `patient_id`, `moh_id`, `age`, `sex`, `os_status`, `os_months`, `smoking_history`. `moh_id` is carried through to the `MOH_ID` column of `data_clinical_patient.txt`.

**Sample file** (tab-separated): `sample_id`, `patient_id`, `cancer_type`, `cancer_type_detailed`, `sample_type`, `primary_tumor_site`, `metastatic_tumor_site`, `tumor_purity`

All listed columns are required. Other columns are ignored. `Unknown` and empty values become `NA`. See `assets/patient_file.txt` and `assets/sample_file.txt` for examples.

> The `sample_id` column in the sample file must use the **deanonymized** (real) sample IDs — the same values that appear in the `deanon_sample_id` column of the linking file. Outputs are scoped to the samplesheet: only samples present in the samplesheet appear in `data_clinical_sample.txt`, `data_clinical_patient.txt`, and `case_lists/`, even if the sample file and patient file contain additional entries.

### 2. Run the pipeline

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results/ \
  --patient_file patient_file.txt \
  --sample_file sample_file.txt \
  --linking_file linking_file.txt \
  --mafsmith_data /path/to/mafsmith_home/ \
  --study_id my_study
```

VCF → MAF uses [mafsmith](https://github.com/nf-osi/mafsmith) + fastVEP (GRCh37, PASS records only), run from the `mafsmith-fastvep` container (see [Containers](#containers)).

`--mafsmith_data` is optional. When it is omitted, empty or `null`, the GRCh37 reference bundle (Ensembl 113) is downloaded once into `assets/mafsmith/` and reused by later runs. A pre-staged bundle must contain `GRCh37/reference.fa` and `GRCh37/genes.gff3.gz`. `MAFSMITH` fails if either file is missing or empty.

> A `-stub-run` with `--skip_mafsmith false` stores an **empty** `assets/mafsmith/`. Delete it afterwards, or later real runs will reuse it and fail on the missing reference files.

Skip VCF → MAF conversion if MAFs already exist:

```bash
nextflow run main.nf ... --skip_mafsmith true
```

Filter mode applies QC thresholds using the sample's `analysis_*_export.tsv`. By default it is off and everything is kept:

```bash
nextflow run main.nf ... --filter_tsv_variants true
```

| Data                            | Kept in filter mode                                                                                              | Params (default)                                                |
| ------------------------------- | ---------------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------- |
| Mutations                       | Position overlaps a TSV row with `Depth` ≥ 250 and `VAF` > 0.03 (`chr` prefix ignored)                           | `--mutation_min_depth` (250), `--mutation_min_vaf` (0.03)       |
| CNAs                            | Raw `Copy Number` ≥ 6, so all deletions are dropped                                                              | `--cna_min_copy_number` (6)                                     |
| Fusions                         | Taken from the TSV, even when a `*-star-fusion.final.vcf` exists. `Supporting Reads` summed per gene pair ≥ 1000 | `--sv_min_supporting_reads` (1000)                              |

With `--fusion_metrics <file.tsv>`, fusions come from that one cohort-wide TSV instead of each sample's export TSV or STAR-Fusion VCF. Its `###` lines are skipped, the `#Sample` line is the header, and rows are matched to samples on `#Sample` (the samplesheet `sample_id`). Only `Total Reads` and `ZR-score` are used as metrics; `Total Reads` becomes `Tumor_Variant_Count`. Fusions are per-sample outputs, so samples already in `--outdir` keep their old fusions: delete `{outdir}/samples/*/*_sv.txt` to rebuild them.
With `--fusion_metrics <file.tsv>`, fusions come from that one cohort-wide TSV instead of each sample's export TSV or STAR-Fusion VCF. Its `###` lines are skipped, the `#Sample` line is the header, and rows are matched to samples on `#Sample`: the samplesheet `sample_id`, the sample folder name, or the folder name without its leading `NNNN_` (folder `1000_23BM004967-1-A01-1` matches `23BM004967-1-A01-1`), ignoring case. Gene names come from `GeneA`/`GeneB`, or from `Fusion` (`ROS1--CD74`) when those hold Ensembl placeholders. Only `Junction Read Count` (reads across the fusion breakpoint) and `ZR-score` are used as metrics; `Junction Read Count` becomes `Tumor_Variant_Count` and `ZR-score` the `FusionMetrics.ZR_score` column of `data_sv.txt` (`NA` for fusions from other sources). `meta_sv.txt` then declares `namespaces: FusionMetrics`, so cBioPortal shows the ZR-score in its structural variant tables. `Total Reads` is ignored: it is the sample's library size, the same for every fusion of a sample. Fusions are per-sample outputs, so samples already in `--outdir` keep their old fusions: delete `{outdir}/samples/*/*_sv.txt` to rebuild them.

`data_seg.txt` is not filtered. Filter mode runs per sample, so samples already in `--outdir` are not re-filtered. Use a fresh `--outdir` when you turn it on or change a threshold.

Optional filters on the merged data. They are off by default and, unlike filter mode, also apply to samples already in `--outdir`:

| Filter       | Drops                                                                                                                                                              | Params (default)                                                                                                         |
| ------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------ |
| gnomAD       | Mutations with a gnomAD allele frequency > 0.5 in any ancestry (exomes or genomes), fetched from Genome Nexus (see [gnomAD frequencies](#gnomad-frequencies))      | `--filter_gnomad_af` (false), `--max_gnomad_af` (0.5), `--gnomad_af_columns` (`gnomAD.AF_max`)                           |
| ClinVar      | Mutations whose ClinVar classification is not Pathogenic, Likely pathogenic or VUS, including those not in ClinVar (see [gnomAD frequencies](#gnomad-frequencies)) | `--filter_clinvar` (false), `--clinvar_keep_classes` (`Pathogenic,Likely_pathogenic,Uncertain_significance`)             |
| VAF          | Mutations with VAF (`t_alt_count / t_depth`) > 0.5                                                                                                                 | `--filter_vaf` (false), `--max_vaf` (0.5)                                                                                |
| CNV depth QC | A sample's CNV calls, unless ≥ 95% of its PASS Pisces variants have depth ≥ 250. The sample is left out of the CNA case lists; mutations and SVs are kept          | `--filter_cnv_depth` (false), `--cnv_min_depth` (250), `--cnv_min_depth_fraction` (0.95), `--cnv_depth_pass_only` (true) |

Values equal to a threshold are kept, and a missing gnomAD or VAF value never drops a mutation; a missing ClinVar classification does. The CNV depth QC writes `qc/cnv_depth_qc.tsv`.

Segment data (`data_seg.txt`, `meta_seg.txt`) is off by default. cBioPortal computes _Fraction Genome Altered_ from it, which is not meaningful for an amplicon panel:

```bash
nextflow run main.nf ... --include_seg true [--seg_mean_mode log2|cn]
```

`log2` (default) writes `log2(CN/2)` and keeps CN=2 segments; `cn` writes the raw copy number and drops CN=2. After changing the mode, delete `{outdir}/samples/*/*_seg.txt`. A study already imported with segment data keeps its FGA values until it is deleted and re-imported.

#### gnomAD frequencies

With `--filter_gnomad_af true` or `--filter_clinvar true`, the pipeline looks up each mutation in [Genome Nexus](https://www.genomenexus.org), the service cBioPortal uses for its own gnomAD and ClinVar columns, and adds these columns to `data_mutations.txt`:

- `gnomADe_AF` (gnomAD v2.1.1 exomes), `gnomADg_AF` (genomes) and `gnomAD.AF` (both combined, the value cBioPortal shows)
- `gnomAD.AF_max`: the highest of the total and the eight ancestry frequencies (afr, amr, asj, eas, fin, nfe, oth, sas), exomes and genomes taken apart; `gnomAD.AF_max_pop` says which one, e.g. `exome_nfe`
- `ClinVar.CLNSIG` (e.g. `Pathogenic/Likely_pathogenic`) and `ClinVar.CLNSIGCONF` (the breakdown of a conflicting call, e.g. `Pathogenic(3)|Uncertain_significance(1)`)

The `gnomAD.*` and `ClinVar.*` columns also appear in the cBioPortal mutation table.

- **gnomAD:** a mutation is dropped when `gnomAD.AF_max` is above `--max_gnomad_af`, so one ancestry above the threshold is enough; one gnomAD has never seen is kept. Use `--gnomad_af_columns gnomAD.AF` for the previous combined-total rule.
- **ClinVar:** a mutation is kept when its classification includes one of `--clinvar_keep_classes` (`Pathogenic/Likely_pathogenic` and `Pathogenic, low penetrance` count). A conflicting call is kept when its breakdown lists one of them. A mutation ClinVar has never seen is **dropped**.
- Genome Nexus sometimes answers without the gnomAD data. Variants missing gnomAD or ClinVar data are asked again (up to two more times) before being treated as absent, which can add up to 30 s to the step.

- The step needs **outbound HTTPS** from the node it runs on. If Genome Nexus can't be reached, the run fails rather than keeping every mutation.
- Only the variant (chromosome, positions, alleles) is sent, never sample or patient IDs.
- To use a self-hosted Genome Nexus (GRCh37), set `--genome_nexus_url`.

Boolean flags accept `true`/`false` on the command line. Nextflow 26 passes CLI values as strings, so the pipeline converts them explicitly.

Resume a previous run:

```bash
nextflow run main.nf ... -resume
```

### Containers

| Label      | Param                  | Default                                                      | Definition                                     |
| ---------- | ---------------------- | ------------------------------------------------------------ | ---------------------------------------------- |
| `python`   | `--python_sif`         | local `.sif`                                                 | `containers/python-ampliseq.def`               |
| `mafsmith` | `--mafsmith_container` | `oras://ghcr.io/crchum-citadel/mafsmith-fastvep:0.1.0-0.4.0` | `containers/mafsmith-fastvep_v0.1.0-0.4.0.def` |

The ghcr images are private; log in first with `apptainer registry login`. To build locally instead, see [`containers/README.md`](containers/README.md).

### Testing

```bash
./nf-test test tests/modules --profile=+apptainer   # all module tests (MAFSMITH pulls the mafsmith image)
nextflow run main.nf -profile test,apptainer        # end-to-end on assets/ test data (skip_mafsmith=true)
```

### Incremental runs

When new samples are added to the samplesheet, re-running with the same `--outdir` automatically skips samples already processed. Only new samples go through per-sample steps; the merge and downstream steps re-run over all samples.

```bash
# Initial run
nextflow run main.nf --input samplesheet_v1.csv --outdir results/ ...

# Later — add new samples to the samplesheet and re-run
nextflow run main.nf --input samplesheet_v2.csv --outdir results/ ...
# → existing samples skipped, new samples processed, outputs updated
```

A sample is skipped when all four files exist under `results/samples/{sample_id}/`:
`{sample_id}_sv.txt`, `{sample_id}_cna.txt`, `{sample_id}_seg.txt`, `{sample_id}_mutations.txt`

### Running standalone scripts (without Nextflow)

The full transformation can also be run via:

```bash
# Edit hardcoded paths at the top of the script first
bash bin/run_pipeline.sh
```

Or run individual scripts from within the output directory:

```bash
cd /path/to/output

python3 /path/to/bin/format_tsv_to_sv.py       <analysis_export.tsv> <SAMPLE_ID>
python3 /path/to/bin/format_cna.py       <analysis_export.tsv> <SAMPLE_ID>
python3 /path/to/bin/format_mutations.py data_mutations.txt    <linking_file>
python3 /path/to/bin/format_sv.py        data_sv.txt           <linking_file>
python3 /path/to/bin/format_cna_deanon.py data_cna.txt         <linking_file>
python3 /path/to/bin/vcf_to_seg.py        <cnv.vcf>             <SAMPLE_ID>
python3 /path/to/bin/seg_deanon.py        data_seg.txt          <linking_file>
python3 /path/to/bin/clinical_patients_format.py <patient_file> <linking_file>
python3 /path/to/bin/clinical_sample_format.py   <sample_file> <linking_file>
```
