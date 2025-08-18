# src/plotter.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_volcano(deg_file, output_path="output/volcano_plot.png"):
    df = pd.read_csv(deg_file)

    # Calculate -log10(p-value)
    df["-log10(p_value)"] = -np.log10(df["p_value"])

    # Set thresholds
    p_thresh = 0.05
    fc_thresh = 1

    # Determine significance
    df["Significant"] = "Not Significant"
    df.loc[(df["p_value"] < p_thresh) & (df["log2FoldChange"].abs() >= fc_thresh), "Significant"] = "Significant"

    # Plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df, x="log2FoldChange", y="-log10(p_value)", hue="Significant", palette={"Significant": "red", "Not Significant": "gray"})
    
    plt.axhline(-np.log10(p_thresh), ls="--", color="black", linewidth=0.7)
    plt.axvline(fc_thresh, ls="--", color="blue", linewidth=0.7)
    plt.axvline(-fc_thresh, ls="--", color="blue", linewidth=0.7)

    plt.title("Volcano Plot of Differential Gene Expression")
    plt.xlabel("log2(Fold Change)")
    plt.ylabel("-log10(p-value)")
    plt.legend()
    plt.tight_layout()
    
    plt.savefig(output_path)
    print(f"✅ Volcano plot saved to: {output_path}")


import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def plot_heatmap(expression_file, deg_file, output_path="output/heatmap_top_genes.png", top_n=10):
    # Load full expression matrix
    expr_df = pd.read_csv(expression_file)
    
    # Load DEGs
    deg_df = pd.read_csv(deg_file)
    
    # Select top N DEGs by p-value
    top_genes = deg_df.sort_values("p_value").head(top_n)["Gene_ID"].tolist()
    
    # Filter expression matrix for top genes
    expr_top = expr_df[expr_df["Gene_ID"].isin(top_genes)].set_index("Gene_ID")
    
    # Normalize data (z-score per gene)
    expr_norm = expr_top.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
    
    # Plot
    plt.figure(figsize=(10, 6))
    sns.heatmap(expr_norm, cmap="coolwarm", linewidths=0.5, annot=True)
    plt.title(f"Top {top_n} Differentially Expressed Genes (Z-score)")
    plt.xlabel("Samples")
    plt.ylabel("Gene ID")
    plt.tight_layout()
    
    plt.savefig(output_path)
    print(f"✅ Heatmap saved to: {output_path}")
