import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as sts
import seaborn as sns
import matplotlib.pyplot as plt
import os

def create_manhattan_plot(results_df, chr_num, results_dir=None):
    # Prepare data for Manhattan plot
    results_df['-log10_p_value'] = -np.log10(results_df['p_value'])

    # Find top 5 windows with highest -log10(p-value)
    top5_windows = results_df.nlargest(5, '-log10_p_value')
    
    plt.figure(figsize=(14, 8))
    
    # Plot all points
    plt.scatter(results_df['win_start'], results_df['-log10_p_value'], c='blue', alpha=0.6, s=30)
    
    # Highlight top 5 windows
    plt.scatter(top5_windows['win_start'], top5_windows['-log10_p_value'], 
                c='red', s=80, alpha=0.8, edgecolors='black', linewidth=1)
    
    # Annotate top 5 windows
    for idx, row in top5_windows.iterrows():
        plt.annotate(
            f"{row['window']}\n-log10(p)={row['-log10_p_value']:.2f}",
            xy=(row['win_start'], row['-log10_p_value']),
            xytext=(10, 10), textcoords='offset points',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
            arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0'),
            fontsize=8, ha='left'
        )
    
    plt.title(f'Manhattan Plot for Chromosome {chr_num}\nTop 5 Most Significant Windows Highlighted')
    plt.xlabel('Genomic Position')
    plt.ylabel('-log10(p-value)')
    plt.axhline(y=-np.log10(0.05), color='red', linestyle='--', alpha=0.7, label='p=0.05 threshold')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Print top 5 windows info
    print(f"\nTop 5 most significant windows for chromosome {chr_num}:")
    print("=" * 60)
    for i, (idx, row) in enumerate(top5_windows.iterrows(), 1):
        print(f"{i}. {row['window']}: -log10(p) = {row['-log10_p_value']:.4f}, "
              f"Position: {row['win_start']:,}-{row['win_end']:,}, "
              f"OR: {row['odds_ratio']:.4f}")
    
    # Determine output path
    if results_dir:
        plot_path = os.path.join(results_dir, f'manhattan_plot_chr{chr_num}.png')
    else:
        plot_path = f'manhattan_plot_chr{chr_num}.png'
    
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Manhattan plot saved to: {plot_path}")
    plt.show()
    
    return top5_windows

if __name__ == "__main__":
    pass