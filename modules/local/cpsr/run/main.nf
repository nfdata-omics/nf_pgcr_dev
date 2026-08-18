process CPSR_RUN {
    tag "${meta.patient}.${meta.sample}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://ghcr.io/sigven/pcgr:2.2.5.singularity'
        : 'nf-core/pcgr:2.2.5'}"

    input:
    tuple val(meta), path(vcf), path(tbi)
    path pcgr_dir
    path vep_cache

    output:
    tuple val(meta), path("${prefix}"), emit: cpsr_reports
    tuple val("${task.process}"), val('cpsr'),  eval("cpsr --version | sed 's/cpsr //g'"), topic: versions, emit: versions_cpsr

    when:
    task.ext.when == null || task.ext.when

    script:
    def genome   = task.ext.genome ?: ''
    def args     = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.id}"
    """
    export XDG_CACHE_HOME=/tmp
    export XDG_DATA_HOME=/tmp
    export QUARTO_PRINT_STACK=true

    mkdir -p ${prefix}

    cpsr \\
        --input_vcf ${vcf} \\
        --vep_dir ${vep_cache} \\
        --refdata_dir ${pcgr_dir} \\
        --output_dir ${prefix} \\
        --genome_assembly ${genome} \\
        --sample_id ${prefix} \\
         ${args}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}
    """
}
