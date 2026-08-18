# nf-core/variantprioritization: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/variantprioritization/usage](https://nf-co.re/variantprioritization/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

nf-core/variantprioritization processes small variant calls (SNVs and InDels) provided in VCF format, with optional copy number alteration (CNA) segments. It generates annotated reports using PCGR for somatic analyses and CPSR for germline analyses.

It is designed to work seamlessly with nf-core/sarek outputs, but also supports VCF files generated independently of Sarek. Supported variant callers include Mutect2 and Strelka for somatic variants, ASCAT for copy numbers calls and DeepVariant or HaplotypeCaller for germline variants.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline.
Use the `--input` parameter to specify its location.

It must be a comma-separated file with a header row and the following required columns:

- `patient`: Patient / subject identifier (no spaces)
- `status`: `0` = normal / germline, `1` = tumor / somatic
- `sample`: Sample identifier (no spaces)
- `vcf`: Path to a VCF / VCF.GZ file

The PCGR / CPSR implementation runs based on the input: If the samplesheet contains both status `0` and `1`, both reporting tools will be used.

Optionally (somatic analysis only), you can provide:

- `cna`: Path to a CNA segments file (required for somatic rows when `--cna_analysis` is enabled)
- `sex`: the patients sex: `MALE`, `FEAMLE`, `UNKNOWN` (default)
- `tumor_site`: the patiens primary tumor site. Use PCGRs numeric codes (see [PCGR docs](https://sigven.github.io/pcgr/articles/running.html#tumor-site) for details). As Tier 1 variants require a tumor site match, PCGR will not find Tier 1 variants if you do not specify a tumor site.
- `tumor_ploidy`: estimated tumor ploidy.
- `tumor_purity`: estimated tumor purity.

```bash
--input '[path to samplesheet file]'
```

### Multiple VCFs per sample

If you have multiple VCFs for the same sample (e.g. one per caller), provide one row per VCF. The pipeline will keep caller-specific files separate during preprocessing and then consolidate/intersect somatic calls per sample for PCGR reporting.

> [!NOTE]
> Intersection of calls is only implemented for the somatic reporting. For germline calls please provide unique sample ids.

### CNA files

If you run with `--cna_analysis`, CNA segments are required for somatic (`status=1`) rows.
The CNA file is expected to be ASCAT-like and will be reformatted into the allele-specific schema required by PCGR.

### Caller detection

The pipeline extracts the variant caller from the VCF header. Make sure to use VCFs that hold this information.

### Example samplesheet

The workflow accepts as input a `samplesheet.csv` file containing the paths to SNV/InDel VCF files and `ASCAT` copy number aberation files. We have efforted to mimic the [samplesheet specifications of nf-core/sarek](https://github.com/nf-core/sarek/blob/master/docs/usage.md#input-sample-sheet-configurations) for ease of use:

| Column  | Description                                                                                                          |
| ------- | :------------------------------------------------------------------------------------------------------------------- |
| patient | Designates the patient/subject; must be unique for each patient, but one patient can have multiple samples           |
| status  | Normal/tumor (0/1) status of sample                                                                                  |
| sample  | Designates the sample ID; must be unique. A patient may have multiple samples e.g a paired tumor-normal, tumor-only. |
| vcf     | Full path to VCF file(s)                                                                                             |
| cna     | Full path to segment file                                                                                            |

An example of a valid samplesheet is given below:

```csv
patient,status,sample,vcf,cna
HCC1395,1,HCC1395T,HCC1395T_vs_HCC1395N.mutect2.vcf.gz,HCC1395T.segments.txt
HCC1395,1,HCC1395T,HCC1395T_vs_HCC1395N.freebayes.vcf.gz,HCC1395T.segments.txt
HCC1395,1,HCC1395T,HCC1395T_vs_HCC1395N.strelka.somatic_snvs.vcf.gz,HCC1395T.segments.txt
HCC1395,1,HCC1395T,HCC1395T_vs_HCC1395N.strelka.somatic_indels.vcf.gz,HCC1395T.segments.txt
HCC1395,0,HCC1395N,HCC1395N.deepvariant.vcf.gz,
HCC1395,0,HCC1395N,HCC1395N.haplotypecaller.vcf.gz,
HCC1396,1,HCC1396T,HCC1396T_vs_HCC1396N.mutect2.vcf.gz,
HCC1396,1,HCC1396T,HCC1396T_vs_HCC1396N.strelka.somatic_snvs.vcf.gz,
HCC1396,1,HCC1396T,HCC1396T_vs_HCC1396N.strelka.somatic_indels.vcf.gz,
```

An example samplesheet is provided in `assets/samplesheet.csv`.

> copy number aberation files must be present for every sample entry when `--cna_analysis true`.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/variantprioritization \
  -profile <docker/singularity/.../institute> \
  --input ./samplesheet.csv \
  --outdir ./results \
  --genome GRCh38
```

> [!WARNING]
> `-profile conda` is not working with the CPSR and PCGR modules at the moment. We are working on a fix and hope to enable it in a future release.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/variantprioritization/usage) and the [parameter documentation](https://nf-co.re/variantprioritization/parameters).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/variantprioritization
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/variantprioritization releases page](https://github.com/nf-core/variantprioritization/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Variant consolidation

Somatic variants called by multiple tools are reformatted to match `PCGR` specifications making them easily searchable in the HTML ouput.

Tumor sample depth (`TDP`), allele frequency (`TAF`) and allelic depths for the ref and alt (`ADT`) are manually calculated and when applicable, applied to the normal sample (`NDP`, `NAF`, `ADN`):

```console
HCC1395T_vs_HCC1395N.mutect2.vcf.gz
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HCC1395_HCC1395N        HCC1395_HCC1395T
chr1    1212740 .       A       C       .       PASS    AS_SB_TABLE=63,80|49,76;DP=282;ECNT=1;MBQ=20,20;MFRL=151,154;MMQ=60,60;MPOS=30;NALOD=1.94;NLOD=25.89;POPAF=6.00;TLOD=341.76     GT:AD:AF:DP:F1R2:F2R1:FAD:SB 0/0:143,0:0.011:143:36,0:36,0:86,0:63,80,0,0    0/1:0,125:0.988:125:0,28:0,37:0,78:0,0,49,76
```

```console
TDP=125;NDP=143;TAF=0.988;NAF=0.011;ADT=0,125;ADN=143,0;TAL=mutect2
```

```console
HCC1395T_vs_HCC1395N.strelka.somatic_snvs.vcf.gz
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NORMAL  TUMOR
chr1    1212740 .       A       C       .       PASS    DP=271;MQ=60.00;MQ0=0;NT=ref;QSS=790;QSS_NT=3070;ReadPosRankSum=0.00;SGT=AA->AC;SNVSB=0.00;SOMATIC;SomaticEVS=19.73;TQSS=1;TQSS_NT=1    DP:FDP:SDP:SUBDP:AU:CU:GU:TU 145:0:0:0:145,145:0,0:0,0:0,0   126:0:0:0:0,0:126,126:0,0:0,0
```

```console
TDP=126;NDP=145;TAF=1;NAF=0;ADT=0,126;ADN=145,0;TAL=strelka
```

Finally, the maximum values for `TAF`, `TDP`, `NAF`, `NDP`, `ADT`, `ADN` are taken as outputs for the consolidate variant call. In addition, values present in the `ID` and `QUAL` column (i.e not `'.'`) are reported if present in any of the original calls:

```console
1       1212740 .       A       C       3793.8  PASS    NDP=145;NAF=0.011;TDP=126;TAF=1;TAL=mutect2,strelka
```

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
