## 1. Take bam files from oncoanalyser test data

### wget samplesheet:

wget https://raw.githubusercontent.com/nf-core/test-datasets/oncoanalyser/samplesheet/bam_eval.subject_a.wgts.tndna_trna.minimal.stub.csv

### wget bams:

wget \
https://raw.githubusercontent.com/nf-core/test-datasets/oncoanalyser/sample_data/simulated_reads/wgts/markdups_bam/subject_a.normal.dna.bwa-mem2_2.2.1.markdups.bam \
https://raw.githubusercontent.com/nf-core/test-datasets/oncoanalyser/sample_data/simulated_reads/wgts/markdups_bam/subject_a.normal.dna.bwa-mem2_2.2.1.markdups.bam.bai \
https://raw.githubusercontent.com/nf-core/test-datasets/oncoanalyser/sample_data/simulated_reads/wgts/markdups_bam/subject_a.tumor.dna.bwa-mem2_2.2.1.markdups.bam \
https://raw.githubusercontent.com/nf-core/test-datasets/oncoanalyser/sample_data/simulated_reads/wgts/markdups_bam/subject_a.tumor.dna.bwa-mem2_2.2.1.markdups.bam.bai

### wget reference:

wget \
https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/GRCh38_masked_exclusions_alts_hlas.fasta \
https://pub-cf6ba01919994c3cbd354659947f74d8.r2.dev/genomes/GRCh38_hmf/25.1/samtools_index-1.16/GRCh38_masked_exclusions_alts_hlas.fasta.fai

## 2. Transform them into vcf.gz files using bcftools.

bcftools mpileup -Ou -f ../reference/GRCh38_masked_exclusions_alts_hlas.fasta ../bam/subject_a.normal.dna.bwa-mem2_2.2.1.markdups.bam | bcftools call -mv -Oz -o subject_a.normal.dna.bwa-mem2_2.2.1.markdups.vcf.gz

bcftools mpileup -Ou -f ../reference/GRCh38_masked_exclusions_alts_hlas.fasta ../bam/subject_a.tumor.dna.bwa-mem2_2.2.1.markdups.bam | bcftool
s call -mv -Oz -o subject_a.tumor.dna.bwa-mem2_2.2.1.markdups.vcf.gz

## . Use those vcf.gz files into pipeline for testing
