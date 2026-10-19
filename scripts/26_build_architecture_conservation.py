import pandas as pd
import numpy as np
from pathlib import Path

# ============================================================
# SCRIPT 26
# ARCHITECTURE REPRESENTATION vs EVOLUTIONARY CONSERVATION
# ============================================================

print("=" * 70)
print("SCRIPT 26: ARCHITECTURE vs CONSERVATION")
print("=" * 70)

# ------------------------------------------------------------
# 1. FILES
# ------------------------------------------------------------

MASTER = Path(
    "data/processed/master_gene_junction_table.tsv"
)

PROFILES = Path(
    "data/processed/correct_phyloP_profiles.tsv"
)

OUTPUT = Path(
    "data/processed/architecture_conservation.tsv"
)

# ------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------

print("\nLoading master table...")
master = pd.read_csv(MASTER, sep="\t")

print("Master rows:", len(master))

print("\nLoading phyloP profiles...")
profiles = pd.read_csv(PROFILES, sep="\t")

print("Profile rows:", len(profiles))

# ------------------------------------------------------------
# 3. CHECK REQUIRED COLUMNS
# ------------------------------------------------------------

master_required = [
    "gene_id",
    "gene_name",
    "chromosome",
    "strand",
    "intron_start",
    "intron_end",
    "architecture_support",
    "total_architectures",
    "architecture_fraction",
]

profile_required = [
    "gene_id",
    "gene_name",
    "strand",
    "intron_start",
    "intron_end",
    "donor_profile",
    "acceptor_profile",
]

for col in master_required:
    if col not in master.columns:
        raise ValueError(f"Missing master column: {col}")

for col in profile_required:
    if col not in profiles.columns:
        raise ValueError(f"Missing profile column: {col}")

# ------------------------------------------------------------
# 4. MATCH PROFILE DATA
# ------------------------------------------------------------

print("\nMatching master table to phyloP profiles...")

profiles_small = profiles[
    [
        "gene_id",
        "strand",
        "intron_start",
        "intron_end",
        "donor_profile",
        "acceptor_profile",
    ]
].copy()

merged = master.merge(
    profiles_small,
    on=[
        "gene_id",
        "strand",
        "intron_start",
        "intron_end",
    ],
    how="inner",
    validate="one_to_one",
)

print("Matched rows:", len(merged))

if len(merged) != len(master):
    print(
        "WARNING:",
        len(master) - len(merged),
        "master rows were not matched."
    )

# ------------------------------------------------------------
# 5. PARSE PHYLOP PROFILES
# ------------------------------------------------------------

POSITIONS = list(range(-20, 0)) + list(range(1, 21))

def parse_profile(value):
    values = [float(x) for x in str(value).split(",")]

    if len(values) != 40:
        raise ValueError(
            f"Expected 40 values, found {len(values)}"
        )

    return values

print("\nParsing phyloP profiles...")

merged["donor_values"] = merged["donor_profile"].apply(parse_profile)
merged["acceptor_values"] = merged["acceptor_profile"].apply(parse_profile)

# ------------------------------------------------------------
# 6. POSITION LOOKUP
# ------------------------------------------------------------

position_index = {
    position: i
    for i, position in enumerate(POSITIONS)
}

def value_at(profile, position):
    return profile[position_index[position]]

# ------------------------------------------------------------
# 7. CALCULATE CONSERVATION FEATURES
# ------------------------------------------------------------

print("\nCalculating conservation features...")

# Donor motif
merged["donor_motif_mean"] = merged["donor_values"].apply(
    lambda x: np.mean([
        value_at(x, 1),
        value_at(x, 2)
    ])
)

# Acceptor motif
merged["acceptor_motif_mean"] = merged["acceptor_values"].apply(
    lambda x: np.mean([
        value_at(x, -2),
        value_at(x, -1)
    ])
)

# Donor splice boundary
merged["donor_boundary_mean"] = merged["donor_values"].apply(
    lambda x: np.mean([
        value_at(x, -2),
        value_at(x, -1),
        value_at(x, 1),
        value_at(x, 2)
    ])
)

# Acceptor splice boundary
merged["acceptor_boundary_mean"] = merged["acceptor_values"].apply(
    lambda x: np.mean([
        value_at(x, -2),
        value_at(x, -1),
        value_at(x, 1),
        value_at(x, 2)
    ])
)

# Complete donor ±20 profile
merged["donor_window_mean"] = merged["donor_values"].apply(
    np.mean
)

# Complete acceptor ±20 profile
merged["acceptor_window_mean"] = merged["acceptor_values"].apply(
    np.mean
)

# ------------------------------------------------------------
# 8. ARCHITECTURE CATEGORY
# ------------------------------------------------------------

merged["architecture_category"] = np.where(
    merged["total_architectures"] == 1,
    "single_architecture_gene",
    "multi_architecture_gene"
)

# ------------------------------------------------------------
# 9. REMOVE TEMPORARY PROFILE COLUMNS
# ------------------------------------------------------------

merged = merged.drop(
    columns=[
        "donor_values",
        "acceptor_values",
        "donor_profile",
        "acceptor_profile",
    ]
)

# ------------------------------------------------------------
# 10. SAVE
# ------------------------------------------------------------

OUTPUT.parent.mkdir(
    parents=True,
    exist_ok=True
)

merged.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)

print("\nSaved:")
print(OUTPUT)

# ------------------------------------------------------------
# 11. SUMMARY
# ------------------------------------------------------------

print("\n" + "=" * 70)
print("SUMMARY")
print("=" * 70)

print("Rows:", len(merged))
print("Genes:", merged["gene_id"].nunique())

print(
    "Single-architecture genes:",
    (merged["total_architectures"] == 1).sum()
)

print(
    "Multi-architecture genes:",
    (merged["total_architectures"] > 1).sum()
)

print(
    "\nArchitecture fraction:"
)

print(
    merged["architecture_fraction"]
    .describe()
    .to_string()
)

print("\nConservation summary:")

print(
    merged[
        [
            "donor_motif_mean",
            "acceptor_motif_mean",
            "donor_boundary_mean",
            "acceptor_boundary_mean",
            "donor_window_mean",
            "acceptor_window_mean",
        ]
    ]
    .describe()
    .to_string()
)

print("\nSCRIPT 26 COMPLETE")
print("=" * 70)
