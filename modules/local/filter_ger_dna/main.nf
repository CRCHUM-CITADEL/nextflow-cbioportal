process FILTER_GER_DNA {
	tag "$meta.sample"	
	label 'process_low'

	input:
	    tuple val(meta), path(ger_dna_vcf)

	output:
	    tuple val(meta), path("*.vcf.gz")

	script:
	"""
	zcat $ger_dna_vcf | grep "#" > tmp.${meta.sample}.vcf
	zcat $ger_dna_vcf | grep PASS | grep -v "FILTER=<ID=low_depth" >> tmp.${meta.sample}.vcf
	bgzip -c tmp.${meta.sample}.vcf > ${meta.sample}.vcf.gz
	"""
}
