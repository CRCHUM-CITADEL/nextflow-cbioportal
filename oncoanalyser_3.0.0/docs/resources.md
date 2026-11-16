# Resource kit

All reference files the pipeline needs, in one versioned download. COSMIC fusions
(`cosmic_data`) are licence-gated and not included.

| Kit version | Pipeline release | Download | Size    | md5 (tar.gz)                       |
| ----------- | ---------------- | -------- | ------- | ---------------------------------- |
| 1.1.0       | 4.1.0            | _TBD_    | _TBD_   | _TBD_                              |
| 1.0.0       | 4.0.0            | _TBD_    | 35.8 GB | `9c8b4d02dc3f89c3c0d095dea3e9c68e` |

## Install and run

```bash
md5sum -c nextflow-cbioportal-resources-v1.1.0.tar.gz.md5
tar -xzf nextflow-cbioportal-resources-v1.1.0.tar.gz
nextflow run main.nf -profile apptainer --resources_dir $PWD/nextflow-cbioportal-resources-v1.1.0 ...
```

## Contents of v1.1.0

| Kit path                                             | Param                                              | Content                                   |
| ---------------------------------------------------- | -------------------------------------------------- | ----------------------------------------- |
| `genomic/annotations/biomart_grch38_…_entrez_id.tsv` | `ensembl_annotations`                              | Ensembl 116 genes with Entrez IDs         |
| `genomic/annotations/ChimerKB4.xlsx`                 | `chimer_data`                                      | ChimerDB 4 known fusions                  |
| `genomic/annotations/cancerhotspots_single.json`     | `cancer_hotspots_data`                             | cancerhotspots.org snapshot               |
| `genomic/vep/cache/`                                 | `vep_data`                                         | VEP 116 cache (used by CPSR)              |
| `genomic/pcgr/`                                      | `pcgr_data`                                        | PCGR bundle 20261106                      |
| `genomic/mafsmith/mafsmith_0.1.0/`                   | `mafsmith_data`                                    | mafsmith home, Ensembl 116, fastVEP 0.5.0 |
| `genomic/reference/Homo_sapiens_assembly38.fasta`    | `genome_reference`                                 | GATK hg38 FASTA (SigProfiler)             |
| `genomic/cosmic_mutational_signatures/`              | `{sbs,dbs,id}_signatures`, `{sbs,dbs,id}_metadata` | COSMIC v3.6 signatures + metadata         |
| `clinical/MoH/dictionary/`                           | `mohccn_*_map`                                     | MOHCCN v3.1 mapping tables                |

Notes:

- PCGR 2.4.x only accepts bundle 20261106 and VEP 116, so all three change together.
- Somatic (mafsmith) and germline (CPSR) rows are both annotated on Ensembl 116, and the
  Entrez table comes from the same release. Change all three together.
- mafsmith uses the kit's fastVEP (0.5.0, as `mafsmith fetch` installs it), not the container's 0.4.0.
  Its transcript cache is pre-built by that fastVEP, so a fastVEP change needs a new kit.
- v1.0.0 → v1.1.0 replaced only the four Ensembl-dependent items above; the rest is unchanged.

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

## Building v1.1.0 from v1.0.0

All Ensembl files come from release 116 (`https://ftp.ensembl.org/pub/release-116/`).

```bash
NEW=nextflow-cbioportal-resources-v1.1.0
cp -r nextflow-cbioportal-resources-v1.0.0 $NEW && cd $NEW

# VEP 116 cache: --vep_dir is the parent of homo_sapiens/116_GRCh38/
rm -rf genomic/vep/cache && mkdir -p genomic/vep/cache
wget https://ftp.ensembl.org/pub/release-116/variation/indexed_vep_cache/homo_sapiens_vep_116_GRCh38.tar.gz
tar -xzf homo_sapiens_vep_116_GRCh38.tar.gz -C genomic/vep/cache && rm homo_sapiens_vep_116_GRCh38.tar.gz

# PCGR bundle 20261106: --refdata_dir holds data/ (check data/grch38/.PCGR_BUNDLE_VERSION)
rm -rf genomic/pcgr && mkdir -p genomic/pcgr
wget https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20261106.grch38.tgz
tar -xzf pcgr_ref_data.20261106.grch38.tgz -C genomic/pcgr && rm pcgr_ref_data.20261106.grch38.tgz

# Ensembl→Entrez table (bin/build_ensembl_annotations.py; BioMart times out on this export)
wget https://ftp.ensembl.org/pub/release-116/gff3/homo_sapiens/Homo_sapiens.GRCh38.116.gff3.gz
wget https://ftp.ensembl.org/pub/release-116/tsv/homo_sapiens/Homo_sapiens.GRCh38.116.entrez.tsv.gz
rm genomic/annotations/biomart_grch38_ensembl_113_with_entrez_id.tsv
build_ensembl_annotations.py Homo_sapiens.GRCh38.116.gff3.gz Homo_sapiens.GRCh38.116.entrez.tsv.gz \
    genomic/annotations/biomart_grch38_ensembl_116_with_entrez_id.tsv

# mafsmith home: Ensembl 116 GFF3 (a real file: `mafsmith fetch --gff3` writes an absolute symlink),
# same reference.fa + .fai, and old fastVEP + its caches out.
M=genomic/mafsmith/mafsmith_0.1.0
rm -rf $M/bin $M/GRCh38/genes.gff3 $M/GRCh38/genes.gff3.fastvep.cache
cp Homo_sapiens.GRCh38.116.gff3.gz $M/GRCh38/genes.gff3.gz

# fastVEP as `mafsmith fetch` installs it (its exact cargo command), built on the image's Debian.
# Plain `mafsmith fetch` would also re-download a 1 GB Ensembl FASTA (it looks for reference.fa.gz).
apptainer exec -B $PWD docker://rust:1.98-bookworm \
    cargo install --git https://github.com/Huang-lab/fastVEP --branch master fastvep-cli --root $M
$M/bin/fastvep --version   # 0.5.0; record it in README.md

# Build fastVEP's transcript cache once, so tasks never write into the kit.
printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t1000000\t.\tA\tG\t.\tPASS\t.\n' > t.vcf
apptainer exec -B $PWD mafsmith-fastvep_v0.1.0-0.4.0.sif bash -c \
    "export HOME=\$(mktemp -d) && ln -s $PWD/$M \$HOME/.mafsmith && mafsmith vcf2maf --input-vcf t.vcf --output-maf t.maf"
ls $M/GRCh38   # genes.gff3  genes.gff3.fastvep.cache  genes.gff3.gz  reference.fa  reference.fa.fai
rm -rf t.vcf t.maf Homo_sapiens.GRCh38.116.*
```

Then steps 2–3 above (`VERSION` = `1.1.0`).
