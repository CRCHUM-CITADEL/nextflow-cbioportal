### mafsmith + fastVEP container (VCF → MAF)

Built from `mafsmith-fastvep_v0.1.0-0.4.0.def` (same image as oncoanalyser):

`apptainer build mafsmith-fastvep_v0.1.0-0.4.0.sif mafsmith-fastvep_v0.1.0-0.4.0.def`

Push it to `oras://ghcr.io/crchum-citadel/mafsmith-fastvep:0.1.0-0.4.0` (the default of `params.mafsmith_container`), or point `--mafsmith_container` at a local `.sif`.

### Python container

Built from `python-ampliseq.def` (`params.python_sif`).
