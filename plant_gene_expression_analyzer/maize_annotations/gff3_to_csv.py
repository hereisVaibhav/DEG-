import pandas as pd

def gff3_to_csv(gff_file, csv_file, subset_size=500):
    """
    Convert GFF3 annotation file into a CSV with Gene_ID and Function.
    Keeps only 'gene' features and extracts descriptions if available.
    
    Args:
        gff_file (str): Path to input .gff3 file
        csv_file (str): Path to output CSV file
        subset_size (int): Number of genes to keep (for MSc project, default = 500)
    """
    genes = []
    with open(gff_file, "r") as f:   # since file is not gzipped
        for line in f:
            if line.startswith("#"):  # skip comments
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            feature_type = parts[2]
            attributes = parts[8]

            # Only keep "gene" features
            if feature_type == "gene":
                attr_dict = {}
                for attr in attributes.split(";"):
                    if "=" in attr:
                        key, value = attr.split("=", 1)
                        attr_dict[key] = value

                gene_id = attr_dict.get("ID", "Unknown")
                description = attr_dict.get("description", "Unknown")

                genes.append({"Gene_ID": gene_id, "Function": description})

    # Convert to DataFrame
    df = pd.DataFrame(genes)

    # Take subset for MSc project (instead of full genome)
    df = df.head(subset_size)

    # Save as CSV
    df.to_csv(csv_file, index=False)
    print(f"✅ Saved annotation subset to {csv_file}")


# =============================
# Example usage
# =============================
if __name__ == "__main__":
    gff3_to_csv(
        "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.62.chromosome.1.gff3",  # Corrected path
        "maize_gene_annotations.csv",
        subset_size=500  # adjust to keep more genes
    )
