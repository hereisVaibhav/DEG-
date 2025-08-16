# main.py

from src.analyzer import load_expression_data, compute_differential_expression

if __name__ == "__main__":
    path = "data/expression_matrix.csv"
    df = load_expression_data(path)

    control_cols = ["Control_1", "Control_2"]
    stress_cols = ["Stress_1", "Stress_2"]

    results = compute_differential_expression(df, control_cols, stress_cols)
    print("\nDifferential Expression Results:\n")
    print(results)

    # Save to output folder
    results.to_csv("output/deg_results.csv", index=False)
    print("\nSaved to output/deg_results.csv")

from src.plotter import plot_volcano, plot_heatmap


# Add this at the end of main.py
plot_volcano("output/deg_results.csv")
plot_heatmap("data/expression_matrix.csv", "output/deg_results.csv")




