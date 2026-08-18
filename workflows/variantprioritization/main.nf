/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { REFERENCE_DATA         } from '../../subworkflows/local/reference_data'

include { INPUT_PREPROCESSING    } from '../../subworkflows/local/input_preprocessing'
include { PREPARE_SOMATIC        } from '../../subworkflows/local/prepare_somatic'
include { PCGR_RUN               } from '../../modules/local/pcgr/run'
include { CPSR_RUN               } from '../../modules/local/cpsr/run'


include { getGenomeAttribute     } from '../../subworkflows/local/utils_nfcore_variantprioritization_pipeline'
include { MULTIQC                } from '../../modules/nf-core/multiqc'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../../subworkflows/local/utils_nfcore_variantprioritization_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow VARIANTPRIORITIZATION {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    params
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:
    params.fasta = getGenomeAttribute('fasta')
    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()
    def fasta = params.fasta ? channel.fromPath(params.fasta, checkIfExists: true).map { it -> [[id: it.baseName], it] }.collect() : channel.empty()

    REFERENCE_DATA(
        params.genome,
        params.pcgr_download,
        params.pcgr_bundleversion,
        params.pcgr_database_dir,
        params.vep_cache,
        params.vep_cache_version,
        params.vep_species,
    )

    def ch_pcgr_dir = REFERENCE_DATA.out.pcgr_dir
    def ch_vep_cache = REFERENCE_DATA.out.vep_cache

    //
    // SUBWORKFLOW: Preprocess VCF files
    //
    INPUT_PREPROCESSING(
        ch_samplesheet,
        fasta,
    )

    //
    // SUBWORKFLOW: Format input files for somatic analyses
    //
    PREPARE_SOMATIC(
        INPUT_PREPROCESSING.out.normalised_somatic,
        INPUT_PREPROCESSING.out.ch_cna_files,
    )

    //
    // SUBWORKFLOW: pcgr (somatic)
    //
    PCGR_RUN(
        PREPARE_SOMATIC.out.pcgr_ready_vcf,
        ch_pcgr_dir.collect(),
        ch_vep_cache.collect(),
    )

    //
    // SUBWORKFLOW: cpsr
    //
    CPSR_RUN(
        INPUT_PREPROCESSING.out.normalised_germline,
        ch_pcgr_dir.collect(),
        ch_vep_cache.collect()
    )


    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [process[process.lastIndexOf(':') + 1..-1], "  ${tool}: ${version}"]
        }
        .groupTuple(by: 0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_' + 'variantprioritization_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )


    //
    // MODULE: MultiQC
    //
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    def ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
        def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))


    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'variantprioritization'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
