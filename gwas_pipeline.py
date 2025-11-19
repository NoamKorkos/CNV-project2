import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sts
import seaborn as sns
import re
import sys
import os
from fisher import create_contingency_tables, run_fisher_test
from manhattan_plot import create_manhattan_plot

def load_data(windows_file):
    windows_df = pd.read_csv(windows_file).dropna(subset=['CNVs'])
    #extract chr number from string
    for index, row in windows_df.iterrows():
        entries = row['CNVs'].split(',')[1:]
        # Extract chromosome numbers using regex
        chr_nums = [int(re.search(r'chr(\d+):(\d+)', entry).group(2)) for entry in entries]
        windows_df.at[index, 'CNVs'] = chr_nums

    return windows_df

def merge_by_cnv(windows_df, unique_counts):
    merged_df = windows_df.merge(unique_counts, left_on='CNVs', right_on='CNV_ID', how='left')
    merged_df = merged_df.drop(columns=['CNV_ID'])

    merged_df = merged_df.groupby('win_num').agg({
        'win_start': 'first',
        'win_end': 'first',
        'CNVs': lambda x: list(x),  # or 'count' if you just want the count
        'unique_count_groupA': 'sum',
        'unique_count_not_groupA': 'sum'
    }).reset_index()

    return merged_df


def load_unique_group_counts(unique_counts_file):
    group_counts = pd.read_csv(unique_counts_file, sep='\t')
    group_counts['CNV_ID'] = group_counts['CNV_ID'].astype(int)

    return group_counts


if __name__ == "__main__":
    chr_num = sys.argv[1]
    
    # Create results directory for this chromosome
    results_dir = f'chr{chr_num}_gwas_results'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        print(f"Created directory: {results_dir}")
    
    windows_df = load_data(f'chro_{chr_num}_with_windows.csv')
    windows_df = windows_df.explode('CNVs')

    unique_counts = load_unique_group_counts(f'Unique_Count_Per_CNV_chr{chr_num}.tsv')
    
    merged_df = merge_by_cnv(windows_df, unique_counts)
    
    # Save merged data to results directory
    merged_output_path = os.path.join(results_dir, f'merged_windows_unique_counts_chr{chr_num}.csv')
    merged_df.to_csv(merged_output_path, index=False)
    print(f"Saved merged data to: {merged_output_path}")

    contingency_tables = create_contingency_tables(chr_num, merged_df)
    results_df = run_fisher_test(contingency_tables, merged_df)
    
    # Save Fisher test results to results directory
    fisher_output_path = os.path.join(results_dir, f'fisher_results_chr{chr_num}.csv')
    results_df.to_csv(fisher_output_path, index=False)
    print(f"Saved Fisher test results to: {fisher_output_path}")

    # Create Manhattan plot and save to results directory
    create_manhattan_plot(results_df, chr_num, results_dir)




