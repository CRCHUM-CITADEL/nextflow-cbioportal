### mafsmith + fastVEP container (VCF → MAF)

`mafsmith-fastvep_v0.1.0-0.4.0.def` builds [mafsmith](https://github.com/nf-osi/mafsmith) 0.1.0 and [fastVEP](https://github.com/Huang-lab/fastVEP) 0.4.0 from source in a Rust build stage. It then copies only the two binaries into a `debian:bookworm-slim` image. It is used by `MAFSMITH` and `DOWNLOAD_MAFSMITH` (process label `mafsmith`).

Build locally. cargo compiles both tools, so this takes a while:

```bash
cd containers/
apptainer build mafsmith-fastvep_v0.1.0-0.4.0.sif mafsmith-fastvep_v0.1.0-0.4.0.def
```

Use the local image instead of the registry one:

```bash
nextflow run main.nf ... --mafsmith_container containers/mafsmith-fastvep_v0.1.0-0.4.0.sif
```

Publish it as the default of `params.mafsmith_container`:

```bash
apptainer push mafsmith-fastvep_v0.1.0-0.4.0.sif oras://ghcr.io/crchum-citadel/mafsmith-fastvep:0.1.0-0.4.0
```

The ghcr package is private. Pulling or pushing needs `apptainer registry login --username <github-user> oras://ghcr.io` with a token that has `read:packages` (or `write:packages` to push).

When bumping either tool, rename the file to `mafsmith-fastvep_v<mafsmith>-<fastvep>.def`, push under the matching tag, and update `params.mafsmith_container` in `nextflow.config`.

### Python container

Built from `python-ampliseq.def` (`params.python_sif`):

```bash
apptainer build python-ampliseq.sif python-ampliseq.def
```
