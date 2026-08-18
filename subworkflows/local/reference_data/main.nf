include { PCGR_GETREF              } from '../../../modules/nf-core/pcgr/getref'
include { UNTAR as UNTAR_PCGR_DB   } from '../../../modules/nf-core/untar'
include { UNTAR as UNTAR_VEP_CACHE } from '../../../modules/nf-core/untar'

workflow REFERENCE_DATA {
    take:
    genome
    pcgr_download
    pcgr_bundleversion
    pcgr_database_dir
    vep_cache
    vep_cache_version
    vep_species

    main:
    if (pcgr_download) {
        def pcgr_genome = genome.toLowerCase()
        PCGR_GETREF([[id: 'pcgr_reference'], pcgr_bundleversion, pcgr_genome])
        ch_pcgr_dir = PCGR_GETREF.out.pcgrref.map { _meta, pcgrref -> pcgrref }
    }
    else {
        // This allows to use a tar.gz file which can also be downsampled for CI testing
        if (file("${pcgr_database_dir}").exists() && pcgr_database_dir.endsWith('.tar.gz')) {
            UNTAR_PCGR_DB([[id: 'pcgr_db'], file("${pcgr_database_dir}")])
            ch_pcgr_dir = UNTAR_PCGR_DB.out.untar.map { _meta, pcgr_db_files -> pcgr_db_files }.collect()
        }
        else {
            ch_pcgr_dir = channel.fromPath(pcgr_database_dir, checkIfExists: true).collect()
        }
    }

    // This allows to use a tar.gz file which can also be downsampled for CI testing
    if (file("${vep_cache}").exists() && vep_cache.endsWith('.tar.gz')) {
        UNTAR_VEP_CACHE([[id: 'vep_cache'], file("${vep_cache}")])
        ch_ensemblvep_cache = UNTAR_VEP_CACHE.out.untar.map { _meta, vep_cache_files -> vep_cache_files }.collect()
    }
    else {
        def vep_genome = genome
        def vep_annotation_cache_key = isCloudUrl(vep_cache) ? "${vep_cache_version}_${vep_genome}/" : ""
        def vep_cache_dir = "${vep_annotation_cache_key}${vep_species}/${vep_cache_version}_${vep_genome}"
        def vep_cache_path_full = file("${vep_cache}/${vep_cache_dir}", type: 'dir')
        if (!vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory()) {
            if (vep_cache == "s3://annotation-cache/vep_cache/") {
                error("This path is not available within annotation-cache.\nPlease check https://annotation-cache.github.io/ to create a request for it.")
            }
            else {
                error("Path provided with VEP cache is invalid.\nMake sure there is a directory named ${vep_cache_dir} in ${vep_cache}.")
            }
        }
        ch_ensemblvep_cache = channel.fromPath(file("${vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
    }

    emit:
    pcgr_dir  = ch_pcgr_dir
    vep_cache = ch_ensemblvep_cache
}

// Helper function to check if cache path is from any cloud provider
def isCloudUrl(cache_url) {
    return cache_url.startsWith("s3://") || cache_url.startsWith("gs://") || cache_url.startsWith("az://")
}
