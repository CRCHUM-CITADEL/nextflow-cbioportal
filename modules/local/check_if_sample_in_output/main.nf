process CHECK_IF_SAMPLE_IN_OUTPUT {
	input:
		path(meta)	

	// will return nothing if sample name is not found in group of output directory
	output:
		eval('test -d "${outputDir}/${meta.group}/${meta.sample}"')
	
	script:
	"""
	true
	"""
}

