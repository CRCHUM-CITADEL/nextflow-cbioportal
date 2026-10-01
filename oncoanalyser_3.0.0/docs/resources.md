# Resource kit

All reference files the pipeline needs, in one versioned download. COSMIC fusions
(`cosmic_data`) are licence-gated and not included.

| Kit version | Pipeline release | Download | Size    | md5 (tar.gz)                       |
| ----------- | ---------------- | -------- | ------- | ---------------------------------- |
| 1.0.0       | 4.0.0            | _TBD_    | 35.8 GB | `9c8b4d02dc3f89c3c0d095dea3e9c68e` |

## Install and run

```bash
md5sum -c nextflow-cbioportal-resources-v1.0.0.tar.gz.md5
tar -xzf nextflow-cbioportal-resources-v1.0.0.tar.gz
nextflow run main.nf -profile apptainer --resources_dir $PWD/nextflow-cbioportal-resources-v1.0.0 ...
```

## Contents of v1.0.0

| Kit path                                             | Param                                              | Content                             |
| ---------------------------------------------------- | -------------------------------------------------- | ----------------------------------- |
| `genomic/annotations/biomart_grch38_…_entrez_id.tsv` | `ensembl_annotations`                              | Ensembl 113 BioMart with Entrez IDs |
| `genomic/annotations/ChimerKB4.xlsx`                 | `chimer_data`                                      | ChimerDB 4 known fusions            |
| `genomic/annotations/cancerhotspots_single.json`     | `cancer_hotspots_data`                             | cancerhotspots.org snapshot         |
| `genomic/vep/cache/`                                 | `vep_data`                                         | VEP 113 cache (used by PCGR)        |
| `genomic/pcgr/`                                      | `pcgr_data`                                        | PCGR bundle 20250314                |
| `genomic/mafsmith/mafsmith_0.1.0/`                   | `mafsmith_data`                                    | mafsmith home, fastVEP 0.3.0        |
| `genomic/reference/Homo_sapiens_assembly38.fasta`    | `genome_reference`                                 | GATK hg38 FASTA (SigProfiler)       |
| `genomic/cosmic_mutational_signatures/`              | `{sbs,dbs,id}_signatures`, `{sbs,dbs,id}_metadata` | COSMIC v3.6 signatures + metadata   |
| `clinical/MoH/dictionary/`                           | `mohccn_*_map`                                     | MOHCCN v3.1 mapping tables          |

Notes:

- PCGR 2.2.x only accepts bundle 20250314 and VEP 113, so all three change together.
- mafsmith uses the kit's fastVEP 0.3.0, not the container's 0.4.0.

## How files are found

Each reference param uses, in order: the value you set, then the kit file under
`resources_dir`, then nothing (VEP/PCGR/mafsmith download, hotspots fetched live).
The kit layout is `resourceLayout()` in `subworkflows/local/utils/main.nf`. The run
stops if `resources_dir/VERSION` differs from `params.resources_version`.

## Releasing a new kit version

1. Copy the previous kit to a new `nextflow-cbioportal-resources-vX.Y.Z` directory and make the change.
2. Update `VERSION` and `README.md`, then regenerate the manifest:
   `find . -type f ! -name MANIFEST.md5 | sort | xargs md5sum > MANIFEST.md5`
3. Pack: `tar -czf <dir>.tar.gz <dir> && md5sum <dir>.tar.gz > <dir>.tar.gz.md5`
4. In the pipeline, bump `resources_version`, update `resourceLayout()` if paths
   changed, and add a row to the table above.
