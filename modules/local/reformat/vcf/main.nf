process REFORMAT_VCF {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3f3833b76564f6d94f8aeedf1b7242920feec139034dc77ef454de8af002a0f9/data'
        : 'community.wave.seqera.io/library/bcftools_pysam_pandas_python:c4549d7814cd0bee'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), emit: vcf
    tuple val(meta), path("${prefix}.vcf.gz.tbi"), emit: tbi, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam
    tuple val("${task.process}"), val('python'), eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def tool = meta.tool.startsWith('strelka') ? 'strelka' : meta.tool
    """
    reformat_vcf.py \\
        --tool ${tool} \\
        --input ${vcf} \\
        --output ${prefix}.vcf.gz
    """

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
