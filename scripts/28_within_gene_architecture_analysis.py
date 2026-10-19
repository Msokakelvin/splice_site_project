import pandas as pd
from scipy.stats import spearmanr
from pathlib import Path

# ============================================================
# SCRIPT 28
# WITHIN-GENE ARCHITECTURE vs CONSERVATION
# ============================================================

print("=" * 70)
print("SCRIPT 28: WITHIN-GENE ARCHITECTURE ANALYSIS")
print("=" * 70)

INPUT = Path(
    "data/processed/architecture_conservation.tsv"
)

df = pd.read_csv(INPUT, sep="\t")

# ------------------------------------------------------------
# MULTI-ARCHITECTURE GENES ONLY
# ------------------------------------------------------------

df = df[
    df["total_architectures"] > 1
].copy()

print("\nMulti-architecture dataset")
print("Rows:", len(df))
print("Genes:", df["gene_id"].nunique())

# ------------------------------------------------------------
# VARIABLES
# ------------------------------------------------------------

variables = {
    "Donor motif": "donor_motif_mean",
    "Acceptor motif": "acceptor_motif_mean",
    "Donor boundary": "donor_boundary_mean",
    "Acceptor boundary": "acceptor_boundary_mean",
    "Donor ±20": "donor_window_mean",
    "Acceptor ±20": "acceptor_window_mean",
}

# ------------------------------------------------------------
# WITHIN-GENE CENTERING
# ------------------------------------------------------------

print("\nCalculating within-gene centered variables...")

gene_mean_arch = df.groupby(
    "gene_id"
)["architecture_fraction"].transform("mean")

df["architecture_fraction_within"] = (
    df["architecture_fraction"]
    - gene_mean_arch
)

# ------------------------------------------------------------
# ANALYSIS
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("WITHIN-GENE SPEARMAN CORRELATIONS")
print("=" * 70)

for name, column in variables.items():

    temp = df[
        [
            "gene_id",
            "architecture_fraction_within",
            column
        ]
    ].dropna()

    rho, p = spearmanr(
        temp["architecture_fraction_within"],
        temp[column]
    )

    print(
        f"{name:20s} "
        f"rho = {rho: .6f}   "
        f"p = {p:.3e}   "
        f"n = {len(temp)}"
    )

# ------------------------------------------------------------
# GENES WITH REAL WITHIN-GENE VARIATION
# ------------------------------------------------------------

variation = df.groupby(
    "gene_id"
)["architecture_fraction"].nunique()

genes_with_variation = variation[
    variation > 1
].index

df_var = df[
    df["gene_id"].isin(genes_with_variation)
].copy()

print("\n" + "=" * 70)
print("GENES WITH >1 DISTINCT ARCHITECTURE FRACTION")
print("=" * 70)

print(
    "Genes with variation:",
    len(genes_with_variation)
)

print(
    "Rows:",
    len(df_var)
)

# ------------------------------------------------------------
# SUMMARY OF WITHIN-GENE ARCHITECTURE VARIATION
# ------------------------------------------------------------

gene_range = df.groupby(
    "gene_id"
)["architecture_fraction"].agg(
    ["min", "max"]
)

gene_range["range"] = (
    gene_range["max"] -
    gene_range["min"]
)

print("\nWithin-gene architecture fraction range:")

print(
    gene_range["range"]
    .describe()
    .to_string()
)

print("\nSCRIPT 28 COMPLETE")
print("=" * 70)
