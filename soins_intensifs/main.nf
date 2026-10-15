nextflow.enable.types = true

include { JOIN_EXPRESSION } from './join_expression'
include { OARFISH_EXPRESSION_TO_CBIOPORTAL } from './oarfish_expression_to_cbioportal'

record Sample {
	sample: String
	quant_file: Path
}

params {
	samplesheet : Channel<Sample> 
	gtf_gz: Path
}

workflow {

	main:
	tpm_ch = OARFISH_EXPRESSION_TO_CBIOPORTAL(params.samplesheet, params.gtf_gz)
	tpm_files = tpm_ch.map { r -> r.tpm}.collect()
	expression = JOIN_EXPRESSION(tpm_files)

	publish:
	tpm = tpm_ch
	expression = expression

}

output {
	tpm {
		path { r -> "${r.sample}" }
	}
	expression {
		path './'
	}
}


