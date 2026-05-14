import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

def plot_volcano(deg_file, output_path):
    """
    Generates a volcano plot from DEG results.
    """
    df = pd.read_csv(deg_file)
    if df.empty:
        return

    plt.figure(figsize=(10, 6))
    
    # Calculate -log10(pvalue) for the Y axis
    df['-log10(pvalue)'] = -np.log10(df['pvalue'] + 1e-10)
    
    # Define significance thresholds
    lfc_thresh = 1.0
    pval_thresh = 0.05
    
    # Color genes based on Significance using Adjusted P-value (FDR)
    df['Significance'] = 'Not Significant'
    
    # Use adj_pvalue for better scientific accuracy if it exists
    sig_col = 'adj_pvalue' if 'adj_pvalue' in df.columns else 'pvalue'

    df.loc[(df['log2FoldChange'] > lfc_thresh) & (df[sig_col] < pval_thresh), 'Significance'] = 'Upregulated'
    df.loc[(df['log2FoldChange'] < -lfc_thresh) & (df[sig_col] < pval_thresh), 'Significance'] = 'Downregulated'
    
    palette = {'Not Significant': 'grey', 'Upregulated': 'red', 'Downregulated': 'blue'}
    
    sns.scatterplot(data=df, x='log2FoldChange', y='-log10(pvalue)', hue='Significance', palette=palette, alpha=0.7)
    
    # Add threshold lines
    plt.axhline(-np.log10(pval_thresh), color='black', linestyle='--', linewidth=1)
    plt.axvline(lfc_thresh, color='black', linestyle='--', linewidth=1)
    plt.axvline(-lfc_thresh, color='black', linestyle='--', linewidth=1)
    
    plt.title('Volcano Plot of Differential Expression')
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 P-Value')
    plt.grid(True, linestyle=':', alpha=0.6)
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_heatmap(expression_file, deg_file, output_path):
    """
    Generates a heatmap for top differentially expressed genes.
    """
    expr_df = pd.read_csv(expression_file)
    deg_df = pd.read_csv(deg_file)
    
    if expr_df.empty or deg_df.empty:
        return

    # Select top 10 genes by absolute log2FoldChange
    deg_df['abs_LFC'] = deg_df['log2FoldChange'].abs()
    top_genes = deg_df.sort_values(by='abs_LFC', ascending=False).head(10)['Gene_ID']
    
    # Filter expression data and set Gene_ID (or Function if available) as index
    heatmap_data = expr_df[expr_df['Gene_ID'].isin(top_genes)].copy()
    
    # Use Gene_ID + Function for labels if Function exists
    if 'Function' in heatmap_data.columns:
        heatmap_data['Label'] = heatmap_data['Gene_ID'] + " (" + heatmap_data['Function'] + ")"
    else:
        heatmap_data['Label'] = heatmap_data['Gene_ID']
        
    heatmap_data = heatmap_data.set_index('Label')
    
    # Select only numeric columns automatically
    heatmap_data = heatmap_data.select_dtypes(include=[np.number])
    
    if heatmap_data.empty:
        print(f"Skipping heatmap: No numeric columns found in {expression_file}")
        return

    # Normalize data for better heatmap visualization (z-score)
    # Handle rows with zero variance or NaN values
    def get_zscore(row):
        std = row.std()
        if std == 0 or np.isnan(std):
            return row - row.mean()
        return (row - row.mean()) / std

    heatmap_z = heatmap_data.apply(get_zscore, axis=1)
    
    if heatmap_z.empty:
        print(f"Skipping heatmap: Z-score calculation resulted in empty data for {expression_file}")
        return
    
    plt.figure(figsize=(12, 10))
    sns.heatmap(heatmap_z, annot=True, cmap='RdBu_r', center=0, cbar_kws={'label': 'Z-Score (Expression)'})
    
    plt.title(f'Expression Heatmap (Top Genes) - {os.path.basename(expression_file)}')
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()

def plot_comparison_bars(comparison_df, output_path):
    """
    Generates a bar plot comparing species stress responses.
    """
    if comparison_df.empty:
        return

    plt.figure(figsize=(14, 8))
    sns.set_style("whitegrid")
    
    plot = sns.barplot(
        data=comparison_df, 
        x='Implied Stress', 
        y='abs_log2FoldChange', 
        hue='Species',
        palette='viridis'
    )
    
    plt.title('Comparison of Stress Response (Absolute Log2 Fold Change)', fontsize=16)
    plt.xlabel('Implied Environmental Stress', fontsize=12)
    plt.ylabel('Magnitude of Response (|Log2FC|)', fontsize=12)
    plt.xticks(rotation=45, ha='right')
    plt.legend(title='Species', loc='upper right')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
