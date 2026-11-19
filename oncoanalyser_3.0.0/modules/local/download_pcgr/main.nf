// Bundle 20261106 is the one PCGR 2.4.x accepts; production uses the kit copy.
process DOWNLOAD_PCGR {
    // Release in the path: an older bundle in assets/ must never be reused.
    storeDir "${projectDir}/assets/pcgr/20261106"

    input:
        val _ready   // gate: only run when there is work to do

    output:
        path "pcgr", emit: data_dir, type: 'dir'

    script:
    """
    mkdir -p pcgr
    wget -P pcgr https://insilico.hpc.uio.no/pcgr/pcgr_ref_data.20261106.grch38.tgz
    tar -zxf pcgr/pcgr_ref_data.20261106.grch38.tgz -C pcgr
    rm pcgr/pcgr_ref_data.20261106.grch38.tgz
    """

    stub:
    """
    mkdir -p pcgr
    """
}
