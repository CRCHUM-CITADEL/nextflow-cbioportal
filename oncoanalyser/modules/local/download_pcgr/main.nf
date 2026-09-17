process DOWNLOAD_PCGR {
    storeDir "${projectDir}/assets/"

    input:
        val _ready   // gate: only run when there is work to do

    output:
        path "pcgr", emit: data_dir, type: 'dir'

    script:
    """
    mkdir -p pcgr
    wget -P pcgr https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20240927.grch38.tgz
    tar -zxf pcgr/pcgr_ref_data.20240927.grch38.tgz -C pcgr
    rm pcgr/pcgr_ref_data.20240927.grch38.tgz
    """

    stub:
    """
    mkdir -p pcgr
    """
}
