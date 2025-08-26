# src/plotter.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

def plot_volcano(deg_file, output_path="output/volcano_plot.png", top_n=5):
    """
    Creates a volcano plot of log2FC vs -log10(p-value).
    Annotates top N significant genes with Gene_ID + Function if available.
    """
    df = pd.read_csv(deg_file)

    # Keep function info if available
    func_map = None
    if "Function" in df.columns:
        func_map = dict(zip(df["Gene_ID"], df["Function"]))
        df = df.drop(columns=["Function"])

    # Calculate -log10(p-value)
    df["-log10(p_value)"] = -np.log10(df["p_value"])

    # Thresholds
    p_thresh = 0.05
    fc_thresh = 1

    # Significant vs not significant
    df["Significant"] = "Not Significant"
    df.loc[(df["p_value"] < p_thresh) & (df["log2FoldChange"].abs() >= fc_thresh),
           "Significant"] = "Significant"

    # Plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=df,
                    x="log2FoldChange", y="-log10(p_value)",
                    hue="Significant",
                    palette={"Significant": "red", "Not Significant": "gray"})

    plt.axhline(-np.log10(p_thresh), ls="--", color="black", linewidth=0.7)
    plt.axvline(fc_thresh, ls="--", color="blue", linewidth=0.7)
    plt.axvline(-fc_thresh, ls="--", color="blue", linewidth=0.7)

    plt.title("Volcano Plot of Differential Gene Expression")
    plt.xlabel("log2(Fold Change)")
    plt.ylabel("-log10(p-value)")
    plt.legend()

    # 🔹 Annotate top N genes
    top_genes = df.nsmallest(top_n, "p_value")
    for _, row in top_genes.iterrows():
        label = row["Gene_ID"]
        if func_map:
            label += f" ({func_map.get(row['Gene_ID'], 'Unknown')})"
        plt.text(row["log2FoldChange"], row["-log10(p_value)"] + 0.2,
                 label, fontsize=8, ha='center')

    plt.tight_layout()
    plt.savefig(output_path)
    print(f"✅ Volcano plot saved to: {output_path}")


def plot_heatmap(expression_file, deg_file,
                 output_path="output/heatmap_top_genes.png", top_n=10):
    """
    Heatmap of top N DEGs across samples.
    Row labels show Gene_ID + Function if available.
    """
    expr_df = pd.read_csv(expression_file)
    deg_df = pd.read_csv(deg_file)

    # Keep Function info if available
    func_map = None
    if "Function" in deg_df.columns:
        func_map = dict(zip(deg_df["Gene_ID"], deg_df["Function"]))
        deg_df = deg_df.drop(columns=["Function"])
    if "Function" in expr_df.columns:
        expr_df = expr_df.drop(columns=["Function"])

    # Select top N genes by p-value
    top_genes = deg_df.sort_values("p_value").head(top_n)["Gene_ID"].tolist()

    # Subset expression data
    expr_top = expr_df[expr_df["Gene_ID"].isin(top_genes)].set_index("Gene_ID")

    # Normalize (z-score per gene)
    expr_norm = expr_top.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    # If function map exists, relabel rows with Gene + Function
    if func_map:
        expr_norm.index = [f"{gene} ({func_map.get(gene, 'Unknown')})" for gene in expr_norm.index]

    # Plot heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(expr_norm, cmap="coolwarm", linewidths=0.5, annot=True)
    plt.title(f"Top {top_n} Differentially Expressed Genes (Z-score)")
    plt.xlabel("Samples")
    plt.ylabel("Gene (Function)")
    plt.tight_layout()

    plt.savefig(output_path)
    print(f"✅ Heatmap saved to: {output_path}")
