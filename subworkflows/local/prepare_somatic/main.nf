/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { REFORMAT_VCF    } from '../../../modules/local/reformat/vcf'
include { REFORMAT_CNA    } from '../../../modules/local/reformat/cna'
include { INTERSECT_VCF   } from '../../../modules/local/intersect/vcf'
include { PCGR_PREPAREVCF } from '../../../modules/local/pcgr/preparevcf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PREPARE_SOMATIC {
    take:
    vcf_files
    cna_files

    main:
    pcgr_header = channel.fromPath("${projectDir}/bin/pcgr_header.txt", checkIfExists: true)

    // Reformat input files
    REFORMAT_VCF(vcf_files)
    REFORMAT_CNA(cna_files)

    vcf_ch = REFORMAT_VCF.out.vcf.join(REFORMAT_VCF.out.tbi)
    cna_ch = params.cna_analysis ? REFORMAT_CNA.out.cna : channel.empty()

    // Intersect somatic variants
    // create master TSV file with variant <-> tool mapping
    // Extract VCF and TBI from channel, choose suitable meta info for merging samples (pop meta.tool, meta.status)
    // < [[ meta.patient, meta.sample], all tool vcfs, all tool tbi ]
    per_sample_somatic = vcf_ch.map { meta, vcf, tbi ->
        def var = [:]
        var.patient = meta.patient
        var.status = meta.status
        var.sample = meta.sample
        var.sex = meta.sex
        var.tumor_site = meta.tumor_site
        var.tumor_purity = meta.tumor_purity
        var.tumor_ploidy = meta.tumor_ploidy
        def tool = meta.tool?.startsWith('strelka') ? 'strelka' : meta.tool
        return [var, vcf, tbi, tool]
    }

    per_sample_somatic_vcfs = per_sample_somatic
        .map { var, vcf, tbi, tool ->
            return [var, vcf, tbi, tool]
        }
        .groupTuple()

    per_sample_somatic_with_tools = per_sample_somatic_vcfs.map { meta, vcfs, tbis, tools ->
        [meta + [tools: tools], vcfs, tbis]
    }

    INTERSECT_VCF(per_sample_somatic_with_tools)


    sample_vcfs_keys = INTERSECT_VCF.out.variant_tool_map
        .map { meta, keys ->
            [meta.subMap(['patient', 'status', 'sample', 'sex', 'tumor_site', 'tumor_purity', 'tumor_ploidy']), keys]
        }
        .join(per_sample_somatic_vcfs)
        .map { meta, keys, vcfs, tbis, _tools -> [meta, keys, vcfs, tbis] }

    PCGR_PREPAREVCF(sample_vcfs_keys, pcgr_header.collect())

    pcgr_ready_vcf = PCGR_PREPAREVCF.out.vcf
        .join(cna_ch, remainder: true)
        .map { meta, vcf, tbi, cna ->
            [meta, vcf, tbi, cna ?: []]
        }

    emit:
    pcgr_ready_vcf
}
