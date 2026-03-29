# crchum-citadel/ampliseq-cbioportal

## Introduction

**crchum-citadel/ampliseq-cbioportal** formats ampliseq genomic data (VCFs + structural variant TSVs) into cBioPortal-compatible files: mutations (MAF), copy number alterations, structural variants, and clinical data.

Requires Nextflow >= 25.04.0.

## Usage

### 1. Prepare input files

**Samplesheet** (`samplesheet.csv`):
```csv
group,subject_id,sample_id,folder_location
cohort_A,PATIENT_001,SAMPLE_001,path/to/samples/SAMPLE_001
```

Each sample folder must contain:
- `analysis_*_export.tsv` — structural variant / CNA export
- `*-basespace-pisces.final.vcf.gz` — compressed VCF

**Linking file** (`linking_file.txt`, tab-separated) — maps anonymized → real IDs:
```
sample_id	deanon_sample_id
SAMPLE_001	PATIENT_001
```

**Patient file** (tab-separated): `patient_id`, `age`, `sex`, `os_status`, `os_months`, `smoking_history`

**Sample file** (tab-separated): `num_id`, `sample_id`, `patient_id`, `cancer_type`, `cancer_type_detailed`, `sample_type`, `tumor_site`, `tumor_purity`

### 2. Run the pipeline

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --outdir results/ \
  --patient_file patient_file.txt \
  --sample_file sample_file.txt \
  --linking_file linking_file.txt \
  --vcf2maf_container community.wave.seqera.io/library/vcf2maf_ensembl-vep:... \
  --vep_data /path/to/vep_data/ \
  --study_id my_study
```

Skip VCF → MAF conversion if MAFs already exist:
```bash
nextflow run main.nf ... --skip_vcf2maf true
```

Pass all mutations through without TSV-coordinate filtering:
```bash
nextflow run main.nf ... --filter_tsv_variants false
```

Resume a previous run:
```bash
nextflow run main.nf ... -resume
```

### Running standalone scripts (without Nextflow)

The full transformation can also be run via:

```bash
# Edit hardcoded paths at the top of the script first
bash bin/run_pipeline.sh
```

Or run individual scripts from within the output directory:

```bash
cd /path/to/output

python3 /path/to/bin/format_tsv.py       <analysis_export.tsv> <SAMPLE_ID>
python3 /path/to/bin/format_cna.py       <analysis_export.tsv> <SAMPLE_ID>
python3 /path/to/bin/format_mutations.py data_mutations.txt    <linking_file>
python3 /path/to/bin/format_sv.py        data_sv.txt           <linking_file>
python3 /path/to/bin/format_cna_deanon.py data_cna.txt         <linking_file>
python3 /path/to/bin/clinical_patients_format.py <patient_file>
python3 /path/to/bin/clinical_sample_format.py   <sample_file>
```

## Credits

crchum-citadel/ampliseq-cbioportal was originally written by Justin.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use crchum-citadel/ampliseq-cbioportal for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->



This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
