import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import os

class CreateDeltaTable:
    def __init__(self, windows_gwas_path, individual_cnv_gwas_path, windows_cnv_join_path):
        self.windows_gwas_path = windows_gwas_path
        self.individual_cnv_gwas_path = individual_cnv_gwas_path
        self.windows_cnv_join_path = windows_cnv_join_path

    def create_delta_table(self):
        windows_gwas_res = pd.read_csv(self.windows_gwas_path) # Roe's results on cnv's. 
        individual_cnv_gwas_res = pd.read_csv(self.individual_cnv_gwas_path) # pvalues per individual cnv
        windows_cnv_join = pd.read_csv(self.windows_cnv_join_path) # which cnv's are in which windows

        #create new table to map cnv to window based off of windows_cnv_join.
        # Explode the CNV list for each window to create individual CNV-to-window mappings
        # assumes each cnv appears in only one window.
        cnv_to_window = windows_cnv_join[['win_num','win_start', 'CNVs']].copy()
        
        # Parse the string representation of lists into actual lists
        import ast
        cnv_to_window['CNVs'] = cnv_to_window['CNVs'].apply(ast.literal_eval)
        cnv_to_window = cnv_to_window.explode('CNVs')
        cnv_to_window = cnv_to_window.rename(columns={'CNVs': 'cnv_id'})
        
        # Convert cnv_id to int to match the format in chr11.csv
        cnv_to_window['cnv_id'] = cnv_to_window['cnv_id'].astype(int)

        # Extract chromosome number from cnv_id in individual_cnv_gwas_res
        individual_cnv_gwas_res['cnv_id'] = individual_cnv_gwas_res['id'].apply(lambda x: int(re.search(r'chr(\d+):(\d+)', x).group(2)))

        # Merge with individual CNV p-values
        cnv_to_window = cnv_to_window.merge(
            individual_cnv_gwas_res[['cnv_id', 'pvalue']], 
            on='cnv_id', 
            how='left'
        )

        # Calculate mean p-value for each window
        cnvs_mean_p = cnv_to_window.groupby(['win_num', 'win_start'])['pvalue'].mean().reset_index()
        cnvs_mean_p = cnvs_mean_p.rename(columns={'pvalue': 'mean_pvalue'})

        # Create windows p-value dataframe
        windows_p = pd.DataFrame({
            'win_id': windows_gwas_res['window'].str.replace('win_', '').astype(int),
            'p_value': windows_gwas_res['p_value']
        })

        # Merge with mean p-values
        windows_p = windows_p.merge(cnvs_mean_p[['win_start', 'win_num', 'mean_pvalue']], left_on='win_id', right_on='win_num', how='left')

        # Calculate delta
        windows_p['delta'] = windows_p['p_value'] - windows_p['mean_pvalue']
        
        # Print debug info
        print(f"Created windows_p with columns: {list(windows_p.columns)}")
        print(f"Windows_p shape: {windows_p.shape}")
        
        # Create the plot
        plt.figure(figsize=(12, 8))
        sns.scatterplot(
            x=windows_p['win_start'],
            y=windows_p['delta'],
            palette='viridis',
            s=90,
            edgecolor='k'
        )

        plt.xlabel('Window Start Position')
        plt.ylabel('Delta (p Windows - p Mean)')
        plt.axhline(y=0.4, color='red', linestyle='--', linewidth=1.5)
        plt.axhline(y=-0.4, color='red', linestyle='--', linewidth=1.5)
        plt.title('Delta between Windows and Mean p-values')

        # Add labels for points with abs(delta) >= 0.4
        for _, row in windows_p[windows_p['delta'].abs() >= 0.4].iterrows():
            plt.text(
                row['win_start'], 
                row['delta'] + 0.02,   # slight vertical offset
                f"{int(row['win_id'])}", 
                ha='center', va='bottom', fontsize=8, color='darkblue'
            )

        plt.show()
        
        return windows_p







    
