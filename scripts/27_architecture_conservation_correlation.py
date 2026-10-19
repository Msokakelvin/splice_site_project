import pandas as pd
from scipy.stats import spearmanr
from pathlib import Path

# ============================================================
# SCRIPT 27
# ARCHITECTURE FRACTION vs CONSERVATION
# EXPLORATORY CORRELATION ANALYSIS
# ============================================================

print("=" * 70)
print("SCRIPT 27: ARCHITECTURE vs CONSERVATION CORRELATION")
print("=" * 70)

INPUT = Path(
    "data/processed/architecture_conservation.tsv"
)

print("\nLoading:")
print(INPUT)

df = pd.read_csv(INPUT, sep="\t")

print("Rows:", len(df))
print("Genes:", df["gene_id"].nunique())

# ------------------------------------------------------------
# ANALYSIS VARIABLES
# ------------------------------------------------------------

analyses = {
    "Donor motif": "donor_motif_mean",
    "Acceptor motif": "acceptor_motif_mean",
    "Donor boundary": "donor_boundary_mean",
    "Acceptor boundary": "acceptor_boundary_mean",
    "Donor ±20": "donor_window_mean",
    "Acceptor ±20": "acceptor_window_mean",
}

# ------------------------------------------------------------
# FUNCTION
# ------------------------------------------------------------

def run_analysis(data, label):

    print("\n" + "=" * 70)
    print(label)
    print("=" * 70)

    print("Rows:", len(data))
    print("Genes:", data["gene_id"].nunique())

    for name, column in analyses.items():

        temp = data[
            [
                "architecture_fraction",
                column
            ]
        ].dropna()

        rho, p = spearmanr(
            temp["architecture_fraction"],
            temp[column]
        )

        print(
            f"{name:20s} "
            f"rho = {rho: .6f}   "
            f"p = {p:.3e}   "
            f"n = {len(temp)}"
        )

# ------------------------------------------------------------
# ANALYSIS 1
# ALL DATA
# ------------------------------------------------------------

run_analysis(
    df,
    "ALL GENE-JUNCTION RECORDS"
)

# ------------------------------------------------------------
# ANALYSIS 2
# MULTI-ARCHITECTURE GENES ONLY
# ------------------------------------------------------------

multi = df[
    df["total_architectures"] > 1
].copy()

run_analysis(
    multi,
    "MULTI-ARCHITECTURE GENES ONLY"
)

# ------------------------------------------------------------
# ANALYSIS 3
# SINGLE-ARCHITECTURE GENES
# ------------------------------------------------------------

single = df[
    df["total_architectures"] == 1
].copy()

print("\n" + "=" * 70)
print("SINGLE-ARCHITECTURE GENES")
print("=" * 70)

print("Rows:", len(single))

print(
    "Unique genes:",
    single["gene_id"].nunique()
)

print(
    "Unique architecture fractions:",
    single["architecture_fraction"].unique()
)

print("\nSCRIPT 27 COMPLETE")
print("=" * 70)
