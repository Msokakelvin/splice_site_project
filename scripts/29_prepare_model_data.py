import pandas as pd
import numpy as np

# ---------------------------------------------------------
# 1. Input and output files
# ---------------------------------------------------------

INPUT = "data/processed/architecture_conservation.tsv"

OUTPUT = "data/processed/model_data.tsv"


# ---------------------------------------------------------
# 2. Load the data
# ---------------------------------------------------------

print("Loading data...")

df = pd.read_csv(INPUT, sep="\t")

print("Original rows:", len(df))


# ---------------------------------------------------------
# 3. Keep only genes with multiple architectures
# ---------------------------------------------------------

df = df[df["total_architectures"] > 1].copy()

print("Rows after multi-architecture filtering:", len(df))


# ---------------------------------------------------------
# 4. Calculate intron length
# ---------------------------------------------------------

df["intron_length"] = (
    df["intron_end"] -
    df["intron_start"] +
    1
)


# ---------------------------------------------------------
# 5. Log-transform intron length
# ---------------------------------------------------------

df["log_intron_length"] = np.log10(
    df["intron_length"]
)


# ---------------------------------------------------------
# 6. Simplify splice class
# ---------------------------------------------------------

df["splice_class_simple"] = np.where(
    df["splice_class"] == "canonical_GT_AG",
    "canonical",
    "noncanonical"
)


# ---------------------------------------------------------
# 7. Keep variables required for modelling
# ---------------------------------------------------------

model_columns = [
    "gene_id",
    "gene_name",
    "junction_id",
    "architecture_fraction",
    "intron_length",
    "log_intron_length",
    "splice_class_simple",
    "donor_motif_mean",
    "acceptor_motif_mean"
]

df_model = df[model_columns].copy()


# ---------------------------------------------------------
# 8. Remove incomplete observations
# ---------------------------------------------------------

df_model = df_model.dropna(
    subset=[
        "gene_id",
        "architecture_fraction",
        "log_intron_length",
        "splice_class_simple",
        "donor_motif_mean",
        "acceptor_motif_mean"
    ]
)


# ---------------------------------------------------------
# 9. Save
# ---------------------------------------------------------

df_model.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# 10. Summary
# ---------------------------------------------------------

print("\nModel dataset created successfully.")
print("Rows:", len(df_model))
print("Unique genes:", df_model["gene_id"].nunique())
print("Unique junctions:", df_model["junction_id"].nunique())

print("\nSplice class:")
print(df_model["splice_class_simple"].value_counts())

print("\nArchitecture fraction:")
print(
    df_model["architecture_fraction"].describe()
)

print("\nIntron length:")
print(
    df_model["intron_length"].describe()
)

print("\nLog10 intron length:")
print(
    df_model["log_intron_length"].describe()
)
