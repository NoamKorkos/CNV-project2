import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import sys

def remove_ids_with_w(cell):
    # Split by comma
    entries = cell.split(',')
    # Keep only entries where the lhs (ID) does NOT contain 'W'
    filtered = [e for e in entries if 'W' not in e.split('=')[0]]
    # Rejoin
    return ','.join(filtered)

def count_genotypes(geno_str):
    geno_str = str(geno_str)
    hetero = geno_str.count("0/1")
    homo = geno_str.count("1/1")
    return pd.Series([hetero, homo])



if __name__ == "__main__":
    chr_num = sys.argv[1]
    path = '/mnt/project/filtered/merged/chr' + chr_num + '.tsv' # from nexus
    df = pd.read_table(path, sep = '\t')

    df["7"] = df["7"].apply(remove_ids_with_w)
    df["7"].str.contains("W")

    cnv_lengths_freq_df = pd.DataFrame({
    "CNV": df['2'],
    "chr num": df['0'],
    "position": df['1'],
    })
    cnv_lengths_freq_df['length'] = df['3'].str.extract(r"SVSIZE=(\d+)").astype(float)

    cnv_lengths_freq_df[["#heterozygote", "#homozygote"]] = df['7'].apply(count_genotypes)
    cnv_lengths_freq_df['quality'] = df['5']

    # Weighted frequency
    cnv_lengths_freq_df["frequency"] = (
        2 * cnv_lengths_freq_df["#homozygote"] +
        cnv_lengths_freq_df["#heterozygote"]
    )

    cnv_lengths_freq_df.to_csv(f'chr{chr_num}_table_for_hadasa_script.csv')
