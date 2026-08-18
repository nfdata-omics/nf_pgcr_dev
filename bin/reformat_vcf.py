#!/usr/bin/env python3

#####################################################################################################################
# Inspired by: @gudeqing
# https://github.com/sigven/pcgr/issues/136#issuecomment-919273152
# Built upon by @BarryDigby to handle Strelka, Freebayes, Mutect2 VCF files.
# Licence: MIT
#####################################################################################################################

from pysam import VariantFile
import os
import argparse
import sys

def mutect2_vaf(record, sample_idx):
    VAF = record.samples[sample_idx]["AF"][0]
    VAF = VAF if VAF is not None else 0
    return VAF


def freebayes_vaf(record, sample_idx):
    AO = record.samples[sample_idx]["AO"][0]
    RO = record.samples[sample_idx]["RO"]
    AO = AO if AO is not None else 0
    RO = RO if RO is not None else 0
    if (AO + RO) == 0:
        VAF = 0
    else:
        VAF = AO / (AO + RO)
    return VAF


# Strelka VCF format:
def strelka_snv_vaf(record, sample_idx):
    ref = str(record.ref + "U")
    alt = str(record.alts[0] + "U")
    tier1RefCounts = record.samples[sample_idx][ref][0]
    tier1AltCounts = record.samples[sample_idx][alt][0]
    tier1RefCounts = tier1RefCounts if tier1RefCounts is not None else 0
    tier1AltCounts = tier1AltCounts if tier1AltCounts is not None else 0
    if (tier1AltCounts + tier1RefCounts) == 0:
        VAF = 0
    else:
        VAF = tier1AltCounts / (tier1AltCounts + tier1RefCounts)
    record.samples[sample_idx]["AD"] = [0, tier1AltCounts]
    return VAF


def strelka_indel_vaf(record, sample_idx):
    tier1RefCounts = record.samples[sample_idx]["TAR"][0]
    tier1AltCounts = record.samples[sample_idx]["TIR"][0]
    tier1RefCounts = tier1RefCounts if tier1RefCounts is not None else 0
    tier1AltCounts = tier1AltCounts if tier1AltCounts is not None else 0
    if (tier1AltCounts + tier1RefCounts) == 0:
        VAF = 0
    else:
        VAF = tier1AltCounts / (tier1AltCounts + tier1RefCounts)
    record.samples[sample_idx]["AD"] = [0, tier1AltCounts]
    return VAF


def strelka_variants_vaf(record, sample_idx):
    AD = record.samples[sample_idx]["AD"]
    ref = AD[0]
    alt = AD[1]
    if (ref + alt) == 0:
        VAF = 0
    else:
        VAF = alt / (ref + alt)
    return VAF


## Strelka TAL = alternative allele, usually the reference
## Strelka TIL = Indel i.e the 'ALT'.
## using tier one [0] as per strelka recommendation.
def strelka_indel_allelic_depth(record, sample_idx):
    tmp_REF = "".join([s.strip() for s in str(record.samples[sample_idx]["TAR"][0])])
    REF = tmp_REF.replace("(", "").replace(")", "")
    tmp_ALT = "".join([s.strip() for s in str(record.samples[sample_idx]["TIR"][0])])
    ALT = tmp_ALT.replace("(", "").replace(")", "")
    return f"{REF},{ALT}"


def strelka_snv_allelic_depth(record, sample_idx):
    ref = str(record.ref + "U")
    alt = str(record.alts[0] + "U")
    tmp_REF = "".join([s.strip() for s in str(record.samples[sample_idx][ref][0])])
    REF = tmp_REF.replace("(", "").replace(")", "")
    tmp_ALT = "".join([s.strip() for s in str(record.samples[sample_idx][alt][0])])
    ALT = tmp_ALT.replace("(", "").replace(")", "")
    return f"{REF},{ALT}"


vcf_formats = {
    "mutect2_vaf": [
        "AD",
        "AF",
        "DP",
        "F1R2",
        "F2R1",
        "FAD",
        "GQ",
        "GT",
        "PGT",
        "PID",
        "PL",
        "PS",
        "SB",
    ],
    "freebayes_vaf": [
        "AD",
        "AO",
        "DP",
        "GL",
        "GQ",
        "GT",
        "MIN_DP",
        "PL",
        "QA",
        "QR",
        "RO",
    ],
    "strelka_snv_vaf": ["AU", "CU", "DP", "FDP", "GU", "SDP", "SUBDP", "TU"],
    "strelka_indel_vaf": [
        "BCN50",
        "DP",
        "DP2",
        "DP50",
        "FDP50",
        "SUBDP50",
        "TAR",
        "TIR",
        "TOR",
    ],
    "strelka_variants_vaf": [
        "AD",
        "ADF",
        "ADR",
        "DP",
        "DPF",
        "DPI",
        "FT",
        "GQ",
        "GQX",
        "GT",
        "MIN_DP",
        "PL",
        "PS",
        "SB",
    ],
}

# --- Flexible Mutect2 detection ---
mutect2_required_formats = {
    "AD",
    "AF",
    "DP",
    "F1R2",
    "F2R1",
    "GQ",
    "GT",
    "PGT",
    "PID",
    "PL",
    "PS",
    "SB",
}

mutect2_optional_formats = {
    "FAD",
}

def detect_vcf_format(header_formats):
    formats_set = set(header_formats)

    # Flexible Mutect2 detection (allow optional FAD)
    if mutect2_required_formats.issubset(formats_set):
        extra = formats_set - mutect2_required_formats - mutect2_optional_formats
        if not extra:
            return "mutect2_vaf"

    # Original strict matching for other callers
    for key, value in vcf_formats.items():
        if set(value) == formats_set:
            return key

    raise ValueError(f"FORMAT combination {list(header_formats)} not recognized.")

def tumor_normal(output_file, tool):
    with VariantFile("tmp_.vcf") as fr:
        header = fr.header
        samples = list(header.samples)
        header.info.add(
            "TDP", number=1, type="Integer", description="Tumor sample depth"
        )
        header.info.add(
            "NDP", number=1, type="Integer", description="Normal sample depth"
        )
        header.info.add("TAF", number=1, type="Float", description="Tumor sample AF")
        header.info.add("NAF", number=1, type="Float", description="Normal sample AF")
        header.info.add(
            "ADT",
            number=".",
            type="String",
            description="Allelic depths for the ref and alt alleles in the order listed (tumor)",
        )
        header.info.add(
            "ADN",
            number=".",
            type="String",
            description="Allelic depths for the ref and alt alleles in the order listed (normal)",
        )
        header.info.add(
            "TAL",
            number=".",
            type="String",
            description="Algorithms that called the somatic mutation",
        )
        samples = list(header.samples)
        formats = list(header.formats)
        fnc_str = detect_vcf_format(formats)
        # end
        header.formats.add(
            "AL",
            number=".",
            type="Integer",
            description="Codes for algorithms that produced the somatic call (1 = freebayes, 2 = mutect2, 3 = strelka)",
        )
        if "strelka" in fnc_str:
            header.formats.add(
                "AD",
                number=2,
                type="Integer",
                description="AD flag for Strelka. Output as tuple so index rule for TAF does not need to be modified.",
            )
        with VariantFile("tmp_1.vcf", "w", header=header) as fw:
            tumor_is_first = 0
            tumor_is_second = 0
            tool_code = 1 if tool == "freebayes" else 2 if tool == "mutect2" else 3
            for record in fr:
                VAF_sample0 = globals()[fnc_str](record, 0)
                VAF_sample1 = globals()[fnc_str](record, 1)
                if VAF_sample0 > VAF_sample1:
                    tumor_is_first += 1
                else:
                    tumor_is_second += 1
                tumor_idx = tumor_is_second >= tumor_is_first
                normal_idx = 1 - tumor_idx
                record.info["TDP"] = record.samples[tumor_idx]["DP"]
                record.info["NDP"] = record.samples[normal_idx]["DP"]
                AF = [VAF_sample0, VAF_sample1]
                record.info["TAF"] = round(AF[tumor_idx], 3)
                record.info["NAF"] = round(AF[normal_idx], 3)
                if fnc_str == "strelka_indel_vaf":
                    record.info["ADT"] = strelka_indel_allelic_depth(record, tumor_idx)
                    record.info["ADN"] = strelka_indel_allelic_depth(record, normal_idx)
                elif fnc_str == "strelka_snv_vaf":
                    record.info["ADT"] = strelka_snv_allelic_depth(record, tumor_idx)
                    record.info["ADN"] = strelka_snv_allelic_depth(record, normal_idx)
                else:
                    tmp = "".join(
                        [s.strip() for s in str(record.samples[tumor_idx]["AD"])]
                    )
                    record.info["ADT"] = tmp.replace("(", "").replace(")", "")
                    tmp = "".join(
                        [s.strip() for s in str(record.samples[normal_idx]["AD"])]
                    )
                    record.info["ADN"] = tmp.replace("(", "").replace(")", "")
                record.info["TAL"] = tool
                record.samples[tumor_idx]["AL"] = tool_code
                record.samples[normal_idx]["AL"] = tool_code
                fw.write(record)
        normal = f"{samples[normal_idx]} NORMAL"
        tumor = f"{samples[tumor_idx]} TUMOR"
        if normal_idx == 0:
            with open("bcftools_reheader.txt", "w") as f:
                f.write(f"{normal}\n{tumor}")
        else:
            with open("bcftools_reheader.txt", "w") as f:
                f.write(f"{tumor}\n{normal}")

    print(f"we guess tumor sample is {samples[tumor_idx]} ")
    out_vcf = "tmp_out.vcf"
    os.system(f"bcftools reheader -s bcftools_reheader.txt tmp_1.vcf > {out_vcf}")
    os.remove("tmp_.vcf")
    os.remove("tmp_1.vcf")
    if ".gz" in output_file:
        os.system(f"bgzip -c {out_vcf} > {output_file}")
        os.remove(out_vcf)
        os.system(f"tabix {output_file}")
    else:
        os.replace(out_vcf, output_file)


def tumor_only(output_file, tool):
    with VariantFile("tmp_.vcf") as fr:
        header = fr.header
        samples = list(header.samples)
        header.info.add(
            "TDP", number=1, type="Integer", description="Tumor sample depth"
        )
        header.info.add("TAF", number=1, type="Float", description="Tumor sample AF")
        header.info.add(
            "ADT",
            number=".",
            type="String",
            description="Allelic depths for the ref and alt alleles in the order listed (tumor)",
        )
        header.info.add(
            "TAL",
            number=".",
            type="String",
            description="Algorithms that called the somatic mutation",
        )
        samples = list(header.samples)
        formats = list(header.formats)
        fnc_str = detect_vcf_format(formats)
        header.formats.add(
            "AL",
            number=".",
            type="Integer",
            description="Codes for algorithms that produced the somatic call (1 = freebayes, 2 = mutect2, 3 = strelka)",
        )
        with VariantFile("tmp_1.vcf", "w", header=header) as fw:
            tumor_idx = 0
            tool_code = 1 if tool == "freebayes" else 2 if tool == "mutect2" else 3
            for record in fr:
                VAF_tumor = globals()[fnc_str](record, 0)
                record.info["TDP"] = record.samples[tumor_idx]["DP"]
                AF = [VAF_tumor]
                record.info["TAF"] = round(AF[tumor_idx], 3)
                tmp = "".join([s.strip() for s in str(record.samples[tumor_idx]["AD"])])
                record.info["ADT"] = tmp.replace("(", "").replace(")", "")
                record.info["TAL"] = tool
                record.samples[tumor_idx]["AL"] = tool_code
                fw.write(record)
        tumor = f"{samples[tumor_idx]} TUMOR"
        with open("bcftools_reheader.txt", "w") as f:
            f.write(f"{tumor}")

    print(f"we guess tumor sample is {samples[tumor_idx]} ")
    out_vcf = "tmp_out.vcf"
    os.system(f"bcftools reheader -s bcftools_reheader.txt tmp_1.vcf > {out_vcf}")
    os.remove("tmp_.vcf")
    os.remove("tmp_1.vcf")
    if ".gz" in output_file:
        os.system(f"bgzip -c {out_vcf} > {output_file}")
        os.remove(out_vcf)
        os.system(f"tabix {output_file}")
    else:
        os.replace(out_vcf, output_file)


def reformat_vcf(vcf_file, output_file, tool):
    os.system(f"bcftools filter -e'FORMAT/DP=\".\"' {vcf_file} -o tmp_.vcf")
    with VariantFile("tmp_.vcf") as fr:
        header = fr.header
        samples = list(header.samples)
        if len(samples) > 1:
            tumor_normal(output_file, tool)
        else:
            tumor_only(output_file, tool)


def main():
    parser = argparse.ArgumentParser(
        description="Reformat VCF files for tumor-normal or tumor-only VCF processing."
    )
    parser.add_argument(
        "-t",
        "--tool",
        required=True,
        choices=["freebayes", "mutect2", "strelka"],
        help="Variant caller used to generate the VCF file.",
    )
    parser.add_argument("-i", "--input", required=True, help="Input VCF file.")
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output file name. If it contains '.gz', bgzip/tabix will be applied.",
    )

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist.")
        sys.exit(1)

    reformat_vcf(args.input, args.output, args.tool)


if __name__ == "__main__":
    main()
