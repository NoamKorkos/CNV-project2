import pandas as pd
import re
import sys
import os

from fisher import create_contingency_tables, run_fisher_test
from manhattan_plot import create_manhattan_plot
from dynamic_windows_creator import DynamicWindowsCreator

def load_data(windows_file):
    windows_df = None
    if 'tsv' in windows_file:
        windows_df = pd.read_csv(windows_file, sep='\t').dropna(subset=['CNV']).reset_index(drop=True)
    else:
        windows_df = pd.read_csv(windows_file).dropna(subset=['CNV']).reset_index(drop=True)
    #extract chr number from string
    for index, row in windows_df.iterrows():
        if ',' in row['CNV']:
            entries = row['CNV'].split(',')[1:]
        else:
            entries = [row['CNV']]

        # Extract chromosome numbers using regex
        chr_nums = [int(re.search(r'chr(\d+):(\d+)', entry).group(2)) for entry in entries]
        if len(chr_nums) == 1:
            windows_df.at[index, 'CNV'] = chr_nums[0]
        else:
            windows_df.at[index, 'CNV'] = chr_nums

    return windows_df

def merge_by_cnv(windows_df, unique_counts):
    merged_df = windows_df.merge(unique_counts, left_on='CNV', right_on='CNV_ID', how='left')
    merged_df = merged_df.drop(columns=['CNV_ID'])

    merged_df = merged_df.groupby('win_num').agg({
        'win_start': 'first',
        'win_end': 'first',
        'CNV': lambda x: list(x),  
        'unique_count_groupA': 'sum',
        'unique_count_not_groupA': 'sum'
    }).reset_index()

    return merged_df


def load_unique_group_counts(unique_counts_file):
    group_counts = pd.read_csv(unique_counts_file, sep='\t').assign(CNV_ID=lambda x: x['CNV_ID'].astype(int))

    return group_counts



if __name__ == "__main__":
    
    ## 1. Call DynamicWIndowCraetor --> pd df
    ## 2. Call Hadasa's fixed window creator
    ## 3. Run GWAS on both sets of windows --> results table, manhattan plot
    ## 4. Call delta 

    chr_num = 11 #sys.argv[1]
    fixed_file_path = '29_10_25/chro_11_with_windows.csv' #sys.argv[2]
    dynamic_file_path = f'29_10_25/chr{chr_num}_full.tsv' #sys.argv[3]

    # Create results directory for this chromosome
    results_dir = f'chr{chr_num}_gwas_results_TEST'
    if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        print(f"Created directory: {results_dir}")

    dynamic_win_df = DynamicWindowsCreator(load_data(dynamic_file_path).assign(CNV=lambda x: x['CNV'].astype(int)))\
                    .adaptive_entropy_windows().explode('CNV').assign(CNV=lambda x: x['CNV'].astype(int))
    dynamic_win_df.to_csv(os.path.join(results_dir, f'dynamic_windows_chr{chr_num}.csv'), index=False)
    fixed_win_df = load_data(fixed_file_path).explode('CNV').assign(CNV=lambda x: x['CNV'].astype(int))

    # Natan's
    unique_counts = load_unique_group_counts(f'29_10_25/Unique_Count_Per_CNV_chr{chr_num}.tsv') 
    

    for df, df_name in [(dynamic_win_df, 'dynamic'), (fixed_win_df, 'fixed')]:
        merged_df = merge_by_cnv(df, unique_counts)

        # Save merged data to results directory
        merged_output_path = os.path.join(results_dir, f'merged_windows_unique_counts_chr{chr_num}_{df_name}.csv')
        merged_df.to_csv(merged_output_path, index=False)
        print(f"Saved merged data to: {merged_output_path}")

        contingency_tables = create_contingency_tables(chr_num, merged_df)
        results_df = run_fisher_test(contingency_tables, merged_df)
        
        # Save Fisher test results to results directory
        fisher_output_path = os.path.join(results_dir, f'fisher_results_chr{chr_num}_{df_name}.csv')
        results_df.to_csv(fisher_output_path, index=False)
        print(f"Saved Fisher test results to: {fisher_output_path}")

        # Create Manhattan plot and save to results directory
        create_manhattan_plot(results_df, chr_num, results_dir, df_name)

        #call Delta
        




