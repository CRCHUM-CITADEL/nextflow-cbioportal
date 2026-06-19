# Output Reference

This document describes every output file produced by the pipeline. All paths are relative to the output directory (`--outdir`).

---

## Genomic mode

### Per-subject files

Intermediate files for each subject are saved in `<outdir>/<group>/<subject>/`:

| File | Description |
|---|---|
| `<sample>_data_cna_hg38.seg` | Copy-number segments (SEG format). One row per genomic segment with chromosome, start, end, number of probes, and segment mean. |
| `<sample>_data_cna_long.txt` | Gene-level discrete copy-number alterations. Columns: Hugo_Symbol, Entrez_Gene_Id, sample, value (-2 to 2 scale). |
| `<sample>.data_sv.txt` | DNA structural variants. Columns include Hugo symbols at both breakpoints, event type (DELETION, DUPLICATION, INVERSION, TRANSLOCATION), SV coordinates, and `DNA_Support=Yes`, `RNA_Support=No`. |
| `<sample>.isofox_fusion.data_sv.txt` | RNA gene fusions from Isofox. Same column schema as DNA SVs but with `Class=FUSION`, `DNA_Support=No`, `RNA_Support=Yes`. Only present for samples with RNA data. |
| `<sample>.tpm.tsv` | Gene expression in TPM (transcripts per million). Columns: Hugo_Symbol, Entrez_Gene_Id, sample TPM value. |
| `<subject>.somatic_rna_germline.maf` | Merged mutation file (MAF format) combining somatic DNA, RNA-supported somatic, and germline variants. Annotated via VEP and PCGR. |
| `<sample>.data_mutational_signatures_contribution_SBS.txt` | SBS (single base substitution) signature contributions fitted against COSMIC v3.6. |
| `<sample>.data_mutational_signatures_counts_SBS.txt` | SBS trinucleotide context counts (96-channel matrix). |
| `<sample>.data_mutational_signatures_contribution_DBS.txt` | DBS (doublet base substitution) signature contributions. |
| `<sample>.data_mutational_signatures_counts_DBS.txt` | DBS dinucleotide context counts (78-channel matrix). |
| `<sample>.data_mutational_signatures_contribution_ID.txt` | ID (indel) signature contributions. |
| `<sample>.data_mutational_signatures_counts_ID.txt` | ID classification counts (83-channel matrix). |

### Per-group aggregated files

All per-subject files are merged into group-level cBioPortal data files in `<outdir>/<group>/`:

#### Data files

| File | cBioPortal type | Description |
|---|---|---|
| `data_cna_hg38.seg` | `COPY_NUMBER_ALTERATION` / `SEG` | Merged copy-number segments for all subjects in the group. |
| `data_cna_long.txt` | `COPY_NUMBER_ALTERATION` / `DISCRETE_LONG` | Merged gene-level discrete CNAs. |
| `data_sv.txt` | `STRUCTURAL_VARIANT` / `SV` | Merged structural variants (DNA SVs + RNA fusions). DNA and RNA rows are distinguished by `DNA_Support`/`RNA_Support` flags. |
| `data_expression.txt` | `MRNA_EXPRESSION` / `CONTINUOUS` | Merged gene expression matrix (TPM). One column per sample, one row per gene. Requires 2+ samples per group. |
| `data_mutations_dna_rna_germline.txt` | `MUTATION_EXTENDED` / `MAF` | Merged mutations in MAF format. Includes somatic DNA, RNA-supported somatic, and germline variants. |
| `data_mutational_signatures_contribution_SBS.txt` | `GENERIC_ASSAY` / `MUTATIONAL_SIGNATURE` | SBS signature contribution matrix. Rows = COSMIC v3.6 SBS signatures, columns = samples. |
| `data_mutational_signatures_counts_SBS.txt` | `GENERIC_ASSAY` / `MUTATIONAL_SIGNATURE` | SBS trinucleotide count matrix. Rows = 96 trinucleotide contexts, columns = samples. |
| `data_mutational_signatures_contribution_DBS.txt` | `GENERIC_ASSAY` / `MUTATIONAL_SIGNATURE` | DBS signature contribution matrix. |
| `data_mutational_signatures_counts_DBS.txt` | `GENERIC_ASSAY` / `MUTATIONAL_SIGNATURE` | DBS dinucleotide count matrix. |
| `data_mutational_signatures_contribution_ID.txt` | `GENERIC_ASSAY` / `MUTATIONAL_SIGNATURE` | ID (indel) signature contribution matrix. |
| `data_mutational_signatures_counts_ID.txt` | `GENERIC_ASSAY` / `MUTATIONAL_SIGNATURE` | ID classification count matrix. |

#### Meta files

Each data file has a corresponding metadata file required by cBioPortal:

| File | Describes |
|---|---|
| `meta_study.txt` | Study-level metadata (name, description, genome build) |
| `meta_cna_hg38.txt` | Copy-number segment data |
| `meta_cna_long.txt` | Gene-level CNA data |
| `meta_sv.txt` | Structural variant data |
| `meta_expression.txt` | Gene expression data |
| `meta_sequenced.txt` | Mutation data |
| `meta_mutational_signatures_contribution_SBS.txt` | SBS signature contributions |
| `meta_mutational_signatures_counts_SBS.txt` | SBS trinucleotide counts |
| `meta_mutational_signatures_contribution_DBS.txt` | DBS signature contributions |
| `meta_mutational_signatures_counts_DBS.txt` | DBS dinucleotide counts |
| `meta_mutational_signatures_contribution_ID.txt` | ID signature contributions |
| `meta_mutational_signatures_counts_ID.txt` | ID classification counts |

#### Case lists

Case lists define which samples belong to each data type, saved in `<outdir>/<group>/case_lists/`:

| File | Description |
|---|---|
| `cases_cnv.txt` | Samples with copy-number data |
| `cases_sequenced.txt` | Samples with mutation data |

#### Utility files

| File | Description |
|---|---|
| `util_linking_file.txt` | TSV mapping `subject_id` to `sample_id` for all subjects in the group. |

#### ML feature tables

Machine-learning-formatted versions of the genomic data, used for downstream analysis:

| File | Description |
|---|---|
| ML-formatted CNV | Gene-level CNAs in ML-compatible format |
| ML-formatted expression | Gene expression in ML-compatible format |
| ML-formatted mutations | Mutations in ML-compatible format |
| ML-formatted SVs | SVs with COSMIC/ChimerDB fusion annotations |

#### Package

| File | Description |
|---|---|
| `<group>.tar.gz` | All cBioPortal data files, meta files, and case lists packaged into a single archive ready for import. |

---

## Clinical mode

Clinical mode outputs are saved in `<outdir>/<group>/`:

| File | Description |
|---|---|
| `data_clinical_patient.txt` | cBioPortal patient attributes. One row per patient with demographics and clinical data. |
| `data_clinical_sample.txt` | cBioPortal sample attributes. One row per sample linked to its patient. |
| `cancer_type.txt` | Cancer type definition file. |
| `meta_clinical_patient.txt` | Metadata for patient attributes. |
| `meta_clinical_sample.txt` | Metadata for sample attributes. |
| `meta_cancer_type.txt` | Metadata for cancer type. |

---

## Pipeline info

Nextflow execution reports are saved in `<outdir>/pipeline_info/`:

| File | Description |
|---|---|
| `execution_report_<timestamp>.html` | Detailed execution report with task-level metrics (CPU, memory, duration). |
| `execution_timeline_<timestamp>.html` | Visual timeline of all tasks. |
| `execution_trace_<timestamp>.txt` | Tab-separated trace of every task (useful for debugging resource usage). |
| `pipeline_dag_<timestamp>.html` | Directed acyclic graph showing the pipeline workflow. |
| `software_versions.yml` | Versions of all tools used in the pipeline run. |
| `params.json` | Parameters used for the pipeline run. |
| `samplesheet.valid.csv` | Validated samplesheet after schema checking. |

---

## Importing into cBioPortal

The pipeline produces a `<group>.tar.gz` archive containing all files needed for cBioPortal import. To import:

1. Extract the archive:

   ```bash
   tar -xzf COHORT1.tar.gz
   ```

2. Validate the study using cBioPortal's validation script:

   ```bash
   python importer/validateData.py -s COHORT1/ -n
   ```

3. Import the study:

   ```bash
   python importer/metaImport.py -s COHORT1/ -u http://localhost:8080 -o
   ```

   Or using the newer importer:

   ```bash
   python importer/cbioportalImporter.py -s COHORT1/ -u http://localhost:8080
   ```

For more details, see the [cBioPortal import documentation](https://docs.cbioportal.org/file-formats/).
