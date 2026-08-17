process INTERSECT_VCF {
    tag "${meta.patient}:${meta.sample}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3f/3f3833b76564f6d94f8aeedf1b7242920feec139034dc77ef454de8af002a0f9/data'
        : 'community.wave.seqera.io/library/bcftools_pysam_pandas_python:c4549d7814cd0bee'}"

    input:
    tuple val(meta), path(vcfs, stageAs: 'inputs/*'), path(tbis, stageAs: 'inputs/*')

    output:
    tuple val(meta), path("${prefix}_keys.txt"), emit: variant_tool_map
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools
    tuple val("${task.process}"), val('pandas'), eval("python -c 'import pandas; print(pandas.__version__)'"), topic: versions, emit: versions_pandas
    tuple val("${task.process}"), val('pysam'), eval("python -c 'import pysam; print(pysam.__version__)'"), topic: versions, emit: versions_pysam
    tuple val("${task.process}"), val('python'), eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    def vcfList = vcfs instanceof List ? vcfs : [vcfs]
    def toolList = meta.tools instanceof List ? meta.tools : [meta.tools]
    def vcfFiles = vcfList.collect { vcfFile -> new File(vcfFile.toString()).name }.join(',')
    def toolNames = toolList.collect { toolName -> toolName.toString() }.join(',')
    """
    isec_vcfs.py \\
        --input_dir inputs \\
        --sample_files ${vcfFiles} \\
        --tool_names ${toolNames} \\
        --output ${prefix}_keys.txt
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_keys.txt
    """
}
