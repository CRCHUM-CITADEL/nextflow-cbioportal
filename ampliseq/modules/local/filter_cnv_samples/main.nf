process FILTER_CNV_SAMPLES {
    tag "${data_file.name}"

    input:
    path(data_file)
    val(failed_ids)   // list of sample IDs that failed CNV depth QC
    val(id_column)    // header name of the sample ID column (e.g. Sample_Id, ID)

    output:
    path("filtered/${data_file.name}")

    script:
    """
    mkdir -p filtered
    # Drop rows whose sample ID failed CNV depth QC (header always kept)
    awk -F'\\t' -v name="${id_column}" -v failed="${failed_ids.join(' ')}" '
        BEGIN  { n = split(failed, ids, " "); for (i = 1; i <= n; i++) skip[ids[i]] = 1 }
        FNR==1 {
            for (i = 1; i <= NF; i++) if (\$i == name) col = i
            if (!col) { print "ERROR: no " name " column in ${data_file.name}" > "/dev/stderr"; exit 1 }
            print; next
        }
        !(\$col in skip)
    ' ${data_file} > "filtered/${data_file.name}"
    """
}
