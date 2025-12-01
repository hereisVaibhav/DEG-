# main.py - UPDATED to handle both AR and Maize 10-gene datasets and comparative plot

from src.analyzer import load_expression_data, compute_differential_expression
from src.plotter import plot_volcano, plot_heatmap, plot_comparison_bars # NEW IMPORT

import pandas as pd
import os

# Static mapping for the 10 synthetic genes' stress types
IMPLIED_STRESS_MAP = {
    "Heat Shock Protein 70": "Heat Stress",
    "Osmotic Stress Regulator (DREB)": "Drought/Osmotic Stress",
    "Salt Tolerance Factor": "Salt Stress",
    "Cold Acclimation Gene (CBF)": "Cold Stress",
    "Jasmonic Acid Defense Signaling": "Biotic/Wounding Stress",
    "Photosynthesis Gene (Rubisco)": "General Stress Suppression",
    "ABA Response Element (AREB)": "Hormone (ABA) Signaling",
    "Iron Storage Protein": "Metal/Toxicity Stress",
    "General Transporter": "Slight Downregulation",
    "Unknown Kinase": "No Significant Change"
}

# Define the input files and the required paths
DATASETS = {
    "AR": {
        "expression_file": "data_matrix/expression_matrix_AR_10_genes_mixed_stress.csv",
        "output_prefix": "AR_10_genes"
    },
    "Maize": {
        "expression_file": "data_matrix/expression_matrix_maize_10_genes_mixed_stress.csv",
        "output_prefix": "Maize_10_genes"
    }
}

# Ensure the output directory exists
os.makedirs("output", exist_ok=True)

def run_analysis_pipeline(species_name, expression_file, output_prefix):
    """Runs the full DGE analysis, saves results, and generates plots for one dataset."""
    print(f"\n--- Starting Analysis for {species_name} ({expression_file}) ---")

    # 1. Load Expression Data
    df = load_expression_data(expression_file)

    # 2. Prepare data for DGE calculation
    functions = None
    string_cols_to_remove = ["Function", "Implied Stress"]
    
    df_for_calc = df.copy()

    # Isolate function info (if present) and remove all non-numeric columns from df_for_calc
    for col in string_cols_to_remove:
        if col in df_for_calc.columns:
            if col == "Function":
                functions = df_for_calc["Function"]
            df_for_calc = df_for_calc.drop(columns=[col])

    control_cols = ["Control_1", "Control_2"]
    stress_cols = ["Stress_1", "Stress_2"]

    # 3. Compute Differential Expression
    results = compute_differential_expression(df_for_calc, control_cols, stress_cols)

    # 4. Add back the function info
    if functions is not None and not results.empty:
        # Merge the functions back using Gene_ID as the common key
        results = results.merge(functions, left_on="Gene_ID", right_index=True)

    # Define output file paths
    deg_output_path = f"output/deg_results_{output_prefix}.csv"
    
    # 5. Save to output folder
    results.to_csv(deg_output_path, index=False)
    print(f"Differential Expression Results for {species_name}:")
    print(results)
    print(f"\n✅ Saved DE results to {deg_output_path}")

    # 6. Generate Plots
    
    # Volcano Plot
    volcano_plot_path = f"output/volcano_plot_{output_prefix}.png"
    plot_volcano(deg_output_path, output_path=volcano_plot_path)
    print(f"✅ Generated Volcano Plot: {volcano_plot_path}")

    # Heatmap
    heatmap_plot_path = f"output/heatmap_top_genes_{output_prefix}.png"
    plot_heatmap(expression_file, deg_output_path, output_path=heatmap_plot_path)
    print(f"✅ Generated Heatmap: {heatmap_plot_path}")


def run_comparison_plot():
    """
    Loads both AR and Maize DE results, creates a combined DataFrame, and plots the comparison.
    """
    all_results = []
    
    print("\n--- Generating Comparison Plot ---")

    for species, info in DATASETS.items():
        deg_file = f"output/deg_results_{info['output_prefix']}.csv"
        
        # Load the saved DE results
        try:
            df = pd.read_csv(deg_file)
        except FileNotFoundError:
            print(f"Skipping comparison: {deg_file} not found. Run individual analysis first.")
            continue

        # Add species and Implied Stress columns
        df['Species'] = species
        # Map Function to Implied Stress
        df['Implied Stress'] = df['Function'].map(IMPLIED_STRESS_MAP)
        # Calculate comparison metric: absolute log2FC
        df['abs_log2FoldChange'] = df['log2FoldChange'].abs()
        
        all_results.append(df)

    if len(all_results) == 2:
        comparison_df = pd.concat(all_results, ignore_index=True)
        comparison_output_path = "output/comparison_bar_plot.png"
        plot_comparison_bars(comparison_df, comparison_output_path)
        print(f"✅ Generated Comparison Bar Plot: {comparison_output_path}")
    else:
        print("Could not generate comparison plot (missing one or both species results).")


if __name__ == "__main__":
    # 1. Run Analysis for AR and Maize individually
    for species, info in DATASETS.items():
        run_analysis_pipeline(species, info["expression_file"], info["output_prefix"])
    
    # 2. Run Comparative Plot
    run_comparison_plot()