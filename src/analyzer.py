# src/analyzer.py

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

def load_expression_data(filepath):
    df = pd.read_csv(filepath)
    df.set_index("Gene_ID", inplace=True)
    return df

def compute_differential_expression(df, control_cols, stress_cols):
    results = []

    for gene in df.index:
        control_vals = df.loc[gene, control_cols].values
        stress_vals = df.loc[gene, stress_cols].values

        # Log2 Fold Change
        mean_control = np.mean(control_vals)
        mean_stress = np.mean(stress_vals)
        log2_fc = np.log2(mean_stress + 1) - np.log2(mean_control + 1)

        # p-value from t-test
        t_stat, p_val = ttest_ind(stress_vals, control_vals)

        results.append({
            "Gene_ID": gene,
            "log2FoldChange": round(log2_fc, 3),
            "p_value": round(p_val, 4)
        })

    result_df = pd.DataFrame(results)
    return result_df
