nextflow.enable.types = true

process OARFISH_EXPRESSION_TO_CBIOPORTAL {

	tag "$sample"
        label 'process_low'

	input:
	record(
		sample: String,
		quant_file: Path
	)
	gtf_gz: Path

	output:
	record(
		sample: sample,
		tpm: file("tpm_${sample}.txt")
	)

	script:
	"""
	gen_oarfish_expression_to_cbioportal.R \\
	--input ${quant_file} \\
	--gtf_gz ${gtf_gz} \\
	--sample_id ${sample} \\
	--output tpm_${sample}.txt
	"""
	
	stub:
	"""
	touch tpm_${sample}.txt
	"""
}
