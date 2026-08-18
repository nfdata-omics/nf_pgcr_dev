# nf-core/variantprioritization: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Reference data](#reference-data) – PCGR and VEP resources when downloaded or provided as archives
- [VCF preprocessing](#vcf-preprocessing) – optional bgzip/tabix, left-normalisation and filtering
- [Variant formatting and intersection](#variant-formatting-and-intersection) – reheader VCFs, merge caller support, and reformat CNAs
- [PCGR](#pcgr) – combined somatic variant annotation and reporting
- [CPSR](#cpsr) – combined germline variant annotation and reporting
- [MultiQC](#multiqc) – aggregate QC and provenance reporting
- [Pipeline information](#pipeline-information) – run metadata and software versions

### Reference data

<details markdown="1">
<summary>Output files</summary>

- `reference/`
  - PCGR reference bundle extracted from `--pcgr_download` or a provided `--pcgr_database_dir` tarball.
  - VEP cache extracted from a provided archive when supplied via `--vep_cache` as a `.tar.gz` bundle.

</details>

Reference resources are only written to the results directory when you request downloads (`--pcgr_download`) or provide archived bundles. When pointing to pre-existing directories, the pipeline reuses them without copying.

### VCF preprocessing

<details markdown="1">
<summary>Output files</summary>

- `htslib/`
  - `*.vcf.gz`: compressed input VCF files
  - `*.tbi`: tabix index files for input VCFs (for example `htslib/HCC1395T_vs_HCC1395N.mutect2.filtered.vcf.gz.tbi`).
- `bcftools/norm/`
  - `{sample}.{caller}.norm.vcf.gz` and `.tbi`: left-aligned, normalized VCFs produced by `bcftools norm` using the reference FASTA.
- `bcftools/filter/`
  - `{sample}.{caller}.norm.filtered.vcf.gz` and `.tbi`: filtered VCFs (using the `bcftools filter` module expression configured in `modules.config`) retaining indexed output.

</details>

Input VCFs are optionally bgzipped and indexed, then normalized and filtered per caller (e.g. Mutect2, Strelka). Output paths retain caller-specific prefixes to make downstream merging transparent.

### Variant formatting and intersection

<details markdown="1">
<summary>Output files</summary>

- `custom/intersect_vcf/`
  - `{sample}_keys.txt`: table mapping each variant to the set of callers that support it (column 5 lists caller names). Used to propagate caller provenance into downstream VCFs.
- `custom/reformat/`
  - `{sample}.{caller}.reformatted.vcf.gz` and `.tbi`: caller VCFs rewritten by `reformat_vcf.py` to add INFO tags such as `TDP`, `NAF`, `TAF`, `ADT`, `ADN`, and caller codes (`TAL`/`AL`). Tumor/normal order is auto-detected and the header is updated accordingly.
  - `{sample}.reformatted.allelic_cna.tsv`: allele-specific CNA table created from CNVkit or ASCAT input via `reformat_cna.py` (`Chromosome`, `Start`, `End`, `nMajor`, `nMinor`).
- `custom/pcgr_ready_vcf/`
  - `{sample}.vcf.gz` and `.tbi`: unified somatic VCF built by combining reformatted caller VCFs and the `{sample}_keys.txt` caller map. This file is ready for PCGR.

</details>

Variants are first intersected across somatic callers to capture support per tool. VCFs are then harmonised (depth/AF fields, caller codes) and merged into a single PCGR-ready VCF per sample. CNAs are reformatted into the allele-specific schema expected by PCGR when CNA analysis is enabled.

### PCGR

<details markdown="1">
<summary>Output files</summary>

- `pcgr/{sample}/`
  - `{sample}.pcgr.{assembly}.html`: interactive PCGR report.
  - `{sample}.pcgr.{assembly}.xlsx`: Excel summary of annotated variants.
  - `{sample}.pcgr.{assembly}.maf`: MAF-formatted somatic variants.
  - `{sample}.pcgr.{assembly}.pass.*`: PASS-filtered TSV/VCF summaries (e.g. `.pass.tsv.gz`, `.pass.vcf.gz` + `.tbi`).
  - `{sample}.pcgr.{assembly}.snv_indel_ann.tsv.gz`, `msigs.tsv.gz`, `tmb.tsv`, `cna_gene*.tsv.gz`, `cna_segment.tsv.gz`: auxiliary annotation tables.
  - `{sample}.pcgr.{assembly}.conf.yaml`: configuration used by PCGR.

</details>

PCGR ingests the unified somatic VCF (and optional allele-specific CNA table) together with the provided reference bundle and VEP cache to produce interactive and machine-readable reports per sample.

### CPSR

<details markdown="1">
<summary>Output files</summary>

- `cpsr/{sample}/`
  - `{sample}.cpsr.{assembly}.html`: interactive CPSR report.
  - `{sample}.cpsr.{assembly}.xlsx`: Excel summary of annotated variants.
  - `{sample}.cpsr.{assembly}.pass.*`: PASS-filtered TSV/VCF summaries (e.g. `.pass.tsv.gz`, `.pass.vcf.gz` + `.tbi`).
  - `{sample}.cpsr.{assembly}.conf.yaml`: configuration used by CPSR.

</details>

CPSR ingests each germline VCF together with the provided reference bundle and VEP cache to produce interactive and machine-readable reports per sample.

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.

### Exemplary Output

#### PCGR

![PCGR](https://raw.githubusercontent.com/sigven/pcgr/master/pcgrr/pkgdown/assets/img/sc2.png)

![PCGR](https://raw.githubusercontent.com/sigven/pcgr/master/pcgrr/pkgdown/assets/img/sc1.png)

![PCGR](https://raw.githubusercontent.com/sigven/pcgr/master/pcgrr/pkgdown/assets/img/sc3.png)

#### CPSR

![CPSR](https://raw.githubusercontent.com/sigven/cpsr/master/pkgdown/assets/img/cpsr_sc.png)
