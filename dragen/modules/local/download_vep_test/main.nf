process DOWNLOAD_VEP_TEST {
    storeDir "${projectDir}/assets"

    output:
        path "vep", emit: cache_dir, type: 'dir'

    script:
    """
    mkdir -p vep
    wget -P vep https://data.cyri.ac/homo_sapiens_vep_112_GRCh38_chr21.tar.gz
    tar -zxf vep/homo_sapiens_vep_112_GRCh38_chr21.tar.gz -C vep
    rm vep/homo_sapiens_vep_112_GRCh38_chr21.tar.gz
    """
}
