# main.py

from src.analyzer import load_expression_data, compute_differential_expression
from src.plotter import plot_volcano, plot_heatmap

import pandas as pd

if __name__ == "__main__":
    path = "data_matrix/expression_matrix_maize_polished.csv"
    df = load_expression_data(path)

    # Remove the Function column from expression data for calculations
    functions = df["Function"]
    df = df.drop(columns=["Function"])

    control_cols = ["Control_1", "Control_2"]
    stress_cols = ["Stress_1", "Stress_2"]

    results = compute_differential_expression(df, control_cols, stress_cols)

    # Add back the function info
    results = results.merge(functions, left_on="Gene_ID", right_index=True)

    print("\nDifferential Expression Results with Functions:\n")
    print(results)

    # Save to output folder
    results.to_csv("output/deg_results_with_functions.csv", index=False)
    print("\n✅ Saved to output/deg_results_with_functions.csv")

    # Plots
    plot_volcano("output/deg_results_with_functions.csv")
    plot_heatmap("data_matrix/expression_matrix_maize_polished.csv",
                 "output/deg_results_with_functions.csv")
