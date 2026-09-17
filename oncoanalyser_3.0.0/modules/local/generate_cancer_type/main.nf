process GENERATE_CANCER_TYPE {
    publishDir { "${params.outdir}/${group}" }, mode: 'copy'

    input:
    val group

    output:
    tuple val(group), path("cancer_type.txt"), path("meta_cancer_type.txt")

    script:
    def type_of_cancer = params.cancer_type ?: params.study_id.toLowerCase().replace('-', '_')
    def name           = params.project_description ?: params.study_id
    """
    printf '%s\\t%s\\t%s\\t%s\\n' '${type_of_cancer}' '${name}' 'Black' 'tissue' > cancer_type.txt
    cat <<EOF > meta_cancer_type.txt
genetic_alteration_type: CANCER_TYPE
datatype: CANCER_TYPE
data_filename: cancer_type.txt
EOF
    """

    stub:
    """
    touch cancer_type.txt
    touch meta_cancer_type.txt
    """
}
