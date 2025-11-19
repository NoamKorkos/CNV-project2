import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sts
import seaborn as sns
import re
import sys

def create_contingency_tables(chr_num, df):
    contingency_tables = {}

    for win_num, group in df.groupby('win_num'):
        # For each window:
        a = group['unique_count_groupA'].sum()         # CNVs in groupA
        b = group['unique_count_not_groupA'].sum()     # CNVs not in groupA

        # Assuming total individuals in each group are known:
        total_groupA = df['unique_count_groupA'].sum()  # adjust if total known separately
        total_not_groupA = df['unique_count_not_groupA'].sum()

        table = pd.DataFrame({
            0: [total_not_groupA - b, total_groupA - a],
            1: [b, a]
        }, index=pd.Index([0, 1], name='groupA'))

        contingency_tables[f'win_{win_num}'] = table

    # Remove any invalid entries if they exist
    if 'win_-1' in contingency_tables:
        contingency_tables.pop('win_-1')

    return contingency_tables


def run_fisher_test(contingency_tables, df):
    # Run Fisher's exact test on each contingency table
    fisher_results = {}
    for win_key, table in contingency_tables.items():
        win_num = int(win_key.split('_')[1])
        OR, pval = sts.fisher_exact(table)

        # Get window start and end positions
        win_start = df[df['win_num'] == win_num]['win_start'].iloc[0]
        win_end = df[df['win_num'] == win_num]['win_end'].iloc[0]
        
        fisher_results[win_key] = {
            'odds_ratio': OR,
            'p_value': pval,
            'win_start': win_start,
            'win_end': win_end
        }

    # Create summary DataFrame
    results_df = pd.DataFrame.from_dict(fisher_results, orient='index')
    results_df.reset_index(inplace=True)
    results_df.rename(columns={'index': 'window'}, inplace=True)

    return results_df

if __name__ == "__main__":
    pass