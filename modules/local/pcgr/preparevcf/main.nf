process PCGR_PREPAREVCF {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/93/935400e4df07528697155f191c225dbf18ac4dc5d7779b1b14f5b974e9237227/data'
        : 'community.wave.seqera.io/library/pysam_vcf2tsvpy_numpy_pandas_python:eb0ee661861e1b56'}"

    input:
    tuple val(meta), path(keys), path(vcfs, stageAs: 'inputs/*'), path(tbis, stageAs: 'inputs/*')
    path pcgr_header

    output:
    tuple val(meta), path("${prefix}.vcf.gz"), path("${prefix}.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('numpy'), eval("python -c 'import numpy; print(numpy.__version__)'"), topic: versions, emit: versions_numpy
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam
    tuple val("${task.process}"), val('python'), eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('vcf2tsvpy'), eval("vcf2tsvpy --version  | cut -d' ' -f2"), topic: versions, emit: versions_vcf2tsvpy

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    pcgr_vcf.py \\
        --sample ${prefix} \\
        --keys ${keys} \\
        --vcf-dir inputs/ \
        --tbi-dir inputs/ \
        --pcgr-header ${pcgr_header}
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    """
}
