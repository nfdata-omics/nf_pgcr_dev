#!/usr/bin/env python3
# Licence: MIT

import numpy as np
import pandas as pd
import argparse
import os
import sys

def reformat_cna(input_file, filename):
    # Read the input file
    df = pd.read_csv(input_file, sep='\t')  # Adjust sep if needed

    #aggiungere uno step in cui leggo header e capisco se è cnvkit o ascat

    if {'sample', 'chr', 'startpos', 'endpos', 'nMajor', 'nMinor'}.issubset(df.columns):
        var = 'ascat'
        df_renamed = df.rename(columns={
            'chr': 'Chromosome',
            'startpos': 'Start',
            'endpos': 'End',
            'nMajor': 'nMajor',
            'nMinor': 'nMinor'
        })
    elif {'chromosome', 'start', 'end', 'cn1', 'cn2'}.issubset(df.columns):
        var = 'cnvkit'

        # Case 1: Fill cn1 = 1, cn2 = cn - 1 when both are missing and cn >= 1
        mask_fill_approx = (
            df['cn1'].isna() & df['cn2'].isna() &
            df['cn'].notna() & np.isfinite(df['cn']) & (df['cn'] >= 1)
        )
        df.loc[mask_fill_approx, 'cn1'] = 1
        df.loc[mask_fill_approx, 'cn2'] = df.loc[mask_fill_approx, 'cn'] - 1

        # Case 2: Fill cn1 = 0, cn2 = 0 when both are missing and cn == 0
        mask_fill_zero = (
            df['cn1'].isna() & df['cn2'].isna() &
            df['cn'].notna() & np.isfinite(df['cn']) & (df['cn'] == 0)
        )
        df.loc[mask_fill_zero, 'cn1'] = 0
        df.loc[mask_fill_zero, 'cn2'] = 0

        df_renamed = df.rename(columns={
            'chromosome': 'Chromosome',
            'start': 'Start',
            'end': 'End',
            'cn1': 'nMajor',
            'cn2': 'nMinor'
        })
    else:
        raise ValueError("Input DataFrame does not match expected column sets for ASCAT or CNVkit.")

    # Select desired columns
    result = df_renamed[['Chromosome', 'Start', 'End', 'nMajor', 'nMinor']]

    # Output filename
    output_file = f"{filename}.tsv"

    # Write to file
    result.to_csv(output_file, sep='\t', index=False)
    print(f"Output written to {output_file}")

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Convert CNVkit .call.cns to PCGR allele-specific format.")
    parser.add_argument("-i", "--input", required=True, help="CNVkit .call.cns input file (tab-delimited).")
    parser.add_argument("-s", "--out", required=True, help="Output filename.")

    args = parser.parse_args()

    # Check if input file exists
    if not os.path.exists(args.input):
        print(f"Error: Input file '{args.input}' does not exist.")
        sys.exit(1)

    # Process CNA file
    reformat_cna(args.input, args.out)

if __name__ == "__main__":
    main()
