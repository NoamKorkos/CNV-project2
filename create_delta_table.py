import pandas as pd
import numpy as np

if __name__ == "__main__":
    windows_gwas_res = pd.read_csv('fisher_results_chr{num}') # fully processed data
    individual_cnv_gwas_res = pd.read_csv('Roei from nexus') # pvalues per individual cnv

    windows_p = pd.DataFrame({
    'win_id': windows_gwas_res['win_id'],
    'p_value': windows_gwas_res['p_value'],
    })

    cnvs_mean_p =...

    windows_p = windows_p.merge(cnvs_mean_p[['win_start', 'win_id']], on='win_id', how='left')
    

    windows_p['delta'] = windows_p['p_value'] - windows_p['mean_pvalue']

    plt.figure(figsize=(8, 5))
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







    
