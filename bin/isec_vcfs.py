#!/usr/bin/env python3
# Licence: MIT

import os
import glob
from pathlib import Path
import pandas as pd
import argparse
import subprocess


def _sort_key_frame(frame):
    frame = frame.copy()
    frame[1] = pd.to_numeric(frame[1], errors="raise")
    return frame.sort_values(by=[0, 1, 2, 3, 4], kind="mergesort").reset_index(
        drop=True
    )


def intersect_variants(input_dir, output_file, sample_files, tool_names):

    input_path = Path(input_dir)
    if not input_path.exists() or not input_path.is_dir():
        raise FileNotFoundError(
            f"Input directory does not exist or is not a directory: {input_dir}"
        )

    if len(sample_files) != len(tool_names):
        raise ValueError(
            "Length mismatch between --sample_files and --tool_names: "
            f"{len(sample_files)} vs {len(tool_names)}"
        )

    if not sample_files:
        raise ValueError(f"No VCF files found in input directory: {input_dir}")

    missing_vcfs = [name for name in sample_files if not (input_path / name).is_file()]
    if missing_vcfs:
        raise FileNotFoundError(
            "Expected staged VCF files are missing in input_dir: "
            + ", ".join(missing_vcfs)
        )

    tbi_files = {
        file.name
        for file in input_path.iterdir()
        if file.is_file() and file.name.endswith(".tbi")
    }
    for vcf_name in sample_files:
        if vcf_name.endswith(".vcf.gz") and f"{vcf_name}.tbi" not in tbi_files:
            raise FileNotFoundError(
                f"Missing index for {vcf_name}: expected {vcf_name}.tbi in {input_dir}"
            )

    tool_map = {idx: tool for idx, tool in enumerate(tool_names)}

    if len(sample_files) > 1:
        for idx, _x in enumerate(sample_files):
            idx = idx + 1  # cant use 0 for -n
            subprocess.run(
                [
                    "bcftools",
                    "isec",
                    "-c",
                    "none",
                    f"-n={idx}",
                    "-p",
                    str(idx),
                    *sample_files,
                ],
                cwd=input_dir,
                check=True,
            )

        pattern = f"{input_dir}/**/sites.txt"
        fn_size = {}
        file_list = glob.glob(pattern, recursive=True)
        for file in file_list:
            file_size = os.stat(file).st_size
            fn_size[file] = file_size
        # remove sites.txt files that are empty
        fn_size = {key: val for key, val in fn_size.items() if val != 0}
        file_list = sorted(fn_size.keys())

        if not file_list:
            raise ValueError(
                "No non-empty sites.txt files produced by bcftools isec. "
                "Check staged VCF/TBI pairs in --input_dir."
            )

        li = []
        for filename in file_list:
            df = pd.read_table(
                filename, sep="\t", header=None, converters={4: str}
            )  # preserve leading zeros
            li.append(df)

        frame = pd.concat(li, axis=0, ignore_index=True)

        # loop over the 4th column containing bytes '0101' etc
        convert_column = []
        for byte_str in frame[4]:
            # print(byte_str)
            # convert 0101 to itemized list
            code = [x for x in byte_str]
            # init list to match 1's in bytestring to corresponding tool name
            grab_index = []
            for idx, val in enumerate(code):
                if val != "0":
                    grab_index.append(int(idx))
            bytes_2_tal = {k: tool_map[k] for k in grab_index if k in tool_map}
            bytes_2_tal = ",".join(
                sorted(bytes_2_tal.values())
            )  # Sort tool names alphabetically
            convert_column.append(bytes_2_tal)

        assert len(convert_column) == len(frame), (
            "bytes to TAL section failed - length of list != length DF"
        )
        frame[4] = convert_column
        # I noticed duplicate rows in the output file during testing. Worrying as I'm not sure how they got there...
        # chr1    3866080 C       T       freebayes
        # chr1    3866080 C       T       freebayes
        frame = frame.drop_duplicates()
        # Sort deterministically across all output columns
        frame = _sort_key_frame(frame)
        frame.to_csv(output_file, sep="\t", index=None, header=None)

    else:
        single_vcf = input_path / sample_files[0]
        result = subprocess.run(
            ["bcftools", "view", str(single_vcf), "-G", "-H"],
            check=True,
            stdout=subprocess.PIPE,
            text=True,
        )
        records = []
        for line in result.stdout.splitlines():
            fields = line.split("\t")
            records.append([fields[0], fields[1], fields[3], fields[4], tool_map[0]])

        frame = pd.DataFrame(records)
        frame = _sort_key_frame(frame)
        frame.to_csv(output_file, sep="\t", index=None, header=None)


def main():
    # Argument parsing using argparse
    parser = argparse.ArgumentParser(
        description="Intersect somatic VCF files for PCGR input."
    )
    parser.add_argument(
        "--input_dir",
        required=True,
        help="Directory containing staged input VCF and TBI files.",
    )
    parser.add_argument(
        "--sample_files",
        required=True,
        help="Comma-separated staged VCF basenames in deterministic order.",
    )
    parser.add_argument(
        "--tool_names",
        required=True,
        help="Comma-separated tool names matching --sample_files order.",
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output key mapping filename."
    )

    args = parser.parse_args()

    # Call intersect_variants function with arguments
    sample_files = [Path(item).name for item in args.sample_files.split(",") if item]
    tool_names = [item for item in args.tool_names.split(",") if item]

    intersect_variants(args.input_dir, args.output, sample_files, tool_names)


if __name__ == "__main__":
    main()
