nextflow.enable.types = true

process JOIN_EXPRESSION {

	label 'process_low'

	input:
	tpm_files: Bag<Path>

	output:
	file('data_expression.txt')

	script:
	"""
	gen_merge_expression_files_to_cbioportal.R \\
	--input_files ${tpm_files.join(',')} \\
	--output_file data_expression.txt \\
	--fill_missing 0 \\
	--strict
	"""

	stub:
	"""
	touch data_expression.txt
	"""
}
