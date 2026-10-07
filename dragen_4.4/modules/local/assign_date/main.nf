process ASSIGN_DATE {

    input:
        tuple val(meta), path(csv)

    output:
        tuple val(meta), path("csv_with_date.csv")

    script:
    """
    cat ${csv}
    awk -F, 'NR==1 {print \$0",date_of_data_extraction"} NR>1 {print \$0",${meta.extraction_date}"}' ${csv} > csv_with_date.csv
    """

}
