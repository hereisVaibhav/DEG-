import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests

def load_expression_data(file_path):
    """
    Loads expression data from a CSV file.
    """
    return pd.read_csv(file_path)

def compute_differential_expression(df, control_cols, stress_cols):
    """
    Computes log2 Fold Change, P-values, and Adjusted P-values (FDR).
    
    Args:
        df (pd.DataFrame): The input expression matrix.
        control_cols (list): List of column names for control samples.
        stress_cols (list): List of column names for stress samples.
        
    Returns:
        pd.DataFrame: A DataFrame with Gene_ID, log2FoldChange, pvalue, and adj_pvalue.
    """
    # Ensure Gene_ID is preserved
    results = df[['Gene_ID']].copy()
    
    # Ensure control and stress columns are numeric
    df[control_cols] = df[control_cols].apply(pd.to_numeric, errors='coerce')
    df[stress_cols] = df[stress_cols].apply(pd.to_numeric, errors='coerce')

    # Calculate means
    mean_control = df[control_cols].mean(axis=1)
    mean_stress = df[stress_cols].mean(axis=1)
    
    # Calculate Log2 Fold Change
    eps = 1e-6
    results['log2FoldChange'] = np.log2((mean_stress + eps) / (mean_control + eps))
    
    # Calculate P-values using Independent T-Test
    def get_pvalue(row_idx):
        ctrl_vals = df.loc[row_idx, control_cols].dropna()
        strs_vals = df.loc[row_idx, stress_cols].dropna()
        
        if len(ctrl_vals) < 2 or len(strs_vals) < 2:
            return 1.0 
            
        _, pval = stats.ttest_ind(ctrl_vals.values.astype(float), strs_vals.values.astype(float))
        return pval

    results['pvalue'] = [get_pvalue(i) for i in df.index]
    results['pvalue'] = results['pvalue'].fillna(1.0)

    # NEW: Calculate Adjusted P-values (FDR) using Benjamini-Hochberg method
    # This is crucial for high-throughput data to control false positives
    _, adj_pvals, _, _ = multipletests(results['pvalue'], method='fdr_bh')
    results['adj_pvalue'] = adj_pvals
    
    return results
