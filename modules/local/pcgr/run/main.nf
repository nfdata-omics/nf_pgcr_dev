process PCGR_RUN {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/sigven/pcgr:2.2.5.singularity'
        : 'nf-core/pcgr:2.2.5'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(cna)
    path pcgr_dir
    path vep_cache

    output:
    tuple val(meta), path("${prefix}"), emit: pcgr_reports
    tuple val("${task.process}"), val('pcgr'),  eval("pcgr --version | sed 's/pcgr //g'"), topic: versions, emit: versions_pcgr

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome           = task.ext.genome      ?: ''
    def args             = task.ext.args        ?: ''
    prefix               = task.ext.prefix      ?: "${meta.id}"
    def cna_cmd          = params.cna_analysis  ? "--input_cna ${cna}"                  : ''
    def opt_sex          = meta.sex             ? "--sex ${meta.sex}"                   : ''
    def opt_tumor_site   = meta.tumor_site      ? "--tumor_site ${meta.tumor_site}"     : ''
    def opt_tumor_purity = meta.tumor_purity    ? "--tumor_purity ${meta.tumor_purity}" : ''
    def opt_tumor_ploidy = meta.tumor_ploidy    ? "--tumor_ploidy ${meta.tumor_ploidy}" : ''
    """
    export XDG_CACHE_HOME=/tmp
    export XDG_DATA_HOME=/tmp
    export QUARTO_PRINT_STACK=true

    mkdir -p ${prefix}

    pcgr \\
        --input_vcf ${vcf} \\
        --vep_dir ${vep_cache} \\
        --refdata_dir ${pcgr_dir} \\
        --output_dir ${prefix} \\
        --genome_assembly ${genome} \\
        --sample_id ${prefix} \\
        --tumor_dp_tag 'TDP' \\
        --tumor_af_tag 'TAF' \\
        --call_conf_tag 'TAL' \\
        ${cna_cmd} \\
        ${opt_sex} \\
        ${opt_tumor_site} \\
        ${opt_tumor_purity} \\
        ${opt_tumor_ploidy} \\
        ${args}
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    """
}
