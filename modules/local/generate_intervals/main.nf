process GENERATE_INTERVALS {
    label 'process_low'

    input:
        path(fai)
        val(chunk_size)

    output:
        path("intervals/*.bed")

    script:
    """
    mkdir intervals
    awk -v chunk=${chunk_size} '{
        chr=\$1; len=\$2
        for(i=0; i<len; i+=chunk) {
            end=(i+chunk>len)?len:i+chunk
            print chr"\\t"i"\\t"end > "intervals/"chr"_"i".bed"
        }
    }' ${fai}
    """
}
