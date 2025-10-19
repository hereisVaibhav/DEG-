# src/analyzer.py

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def load_expression_data(filepath):
    """
    Loads expression data CSV into a DataFrame.
    Keeps 'Function' column if present, sets Gene_ID as index.
    """
    df = pd.read_csv(filepath)
    df.set_index("Gene_ID", inplace=True)
    return df

def compute_differential_expression(df, control_cols, stress_cols):
    """
    Computes log2 fold change and p-values for each gene.
    Works even if 'Function' column is present (ignored during calc).
    """
    results = []

    # If 'Function' exists, separate it out
    functions = None
    if "Function" in df.columns:
        functions = df["Function"]
        df = df.drop(columns=["Function"])

    for gene in df.index:
        control_vals = df.loc[gene, control_cols].values
        stress_vals = df.loc[gene, stress_cols].values

        # Calculate log2 fold change
        mean_control = np.mean(control_vals)
        mean_stress = np.mean(stress_vals)
        log2_fc = np.log2(mean_stress + 1) - np.log2(mean_control + 1)

        # Statistical significance (Welch's t-test)
        t_stat, p_val = ttest_ind(stress_vals, control_vals, equal_var=False)

        results.append({
            "Gene_ID": gene,
            "log2FoldChange": round(log2_fc, 3),
            "p_value": round(p_val, 4)
        })

    result_df = pd.DataFrame(results)

    # Add back functions if available
    if functions is not None:
        result_df = result_df.merge(functions, left_on="Gene_ID", right_index=True)

    return result_df
