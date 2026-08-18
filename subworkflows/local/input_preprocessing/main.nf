/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { BCFTOOLS_FILTER   } from '../../../modules/nf-core/bcftools/filter'
include { BCFTOOLS_NORM     } from '../../../modules/nf-core/bcftools/norm'
include { HTSLIB_BGZIPTABIX } from '../../../modules/nf-core/htslib/bgziptabix/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow INPUT_PREPROCESSING {
    take:
    ch_samplesheet
    fasta

    main:
    def ch_vcf_to_bgziptabix = ch_samplesheet.map { meta, vcf, tbi, _cna -> [meta, vcf, tbi, []] }

    ch_cna_files = ch_samplesheet
        .map { meta, _vcf, _tbi, cna ->
            [meta.subMap(['patient', 'status', 'sample', 'sex', 'tumor_site', 'tumor_purity', 'tumor_ploidy']), cna]
        }
        .distinct()

    HTSLIB_BGZIPTABIX(ch_vcf_to_bgziptabix, 'compress', true, 'vcf')
    def ch_vcf_preprocessed = HTSLIB_BGZIPTABIX.out.output.join(HTSLIB_BGZIPTABIX.out.index)

    BCFTOOLS_NORM(ch_vcf_preprocessed, fasta)
    ch_norm = BCFTOOLS_NORM.out.vcf.join(BCFTOOLS_NORM.out.index)

    BCFTOOLS_FILTER(ch_norm)
    ch_filtered = BCFTOOLS_FILTER.out.vcf.join(BCFTOOLS_FILTER.out.index)

    normalised_germline = ch_filtered.filter { meta, _vcf, _tbi -> meta.status == 'germline' }
    normalised_somatic = ch_filtered.filter { meta, _vcf, _tbi -> meta.status == 'somatic' }

    emit:
    normalised_germline
    normalised_somatic
    ch_cna_files
}
