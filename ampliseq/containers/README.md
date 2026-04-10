### In order to create the vcf2maf container, you need to have wave installed.

`wave --conda vcf2maf=1.6.22=hdfd78af_2 --conda-package ensembl-vep=113.4 --freeze`

### This will create a link to a docker repo like so to use in the module (usable with apptainer): 

community.wave.seqera.io/library/vcf2maf_ensembl-vep:1b486a30e76e2908
