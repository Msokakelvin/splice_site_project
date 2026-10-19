import pandas as pd


# ============================================================
# SCRIPT 16
# Build clean master gene-junction dataset
# ============================================================

FEATURE_FILE = "data/processed/splice_site_features.tsv"
ARCH_FILE = "data/processed/splice_site_architecture_support.tsv"

OUTPUT_FILE = "data/processed/master_gene_junction_table.tsv"


# ------------------------------------------------------------
# 1. Load the feature table
# ------------------------------------------------------------

features = pd.read_csv(
    FEATURE_FILE,
    sep="\t"
)

print("Feature rows:", len(features))


# ------------------------------------------------------------
# 2. Convert intron coordinates to junction coordinates
# ------------------------------------------------------------

def make_junction_id(row):

    if row["strand"] == "+":

        donor = int(row["intron_start"]) - 1
        acceptor = int(row["intron_end"]) + 1

    else:

        donor = int(row["intron_end"]) + 1
        acceptor = int(row["intron_start"]) - 1

    return (
        f'{row["chromosome"]}:'
        f'{donor}-{acceptor}:'
        f'{row["strand"]}'
    )


features["junction_id"] = features.apply(
    make_junction_id,
    axis=1
)


# ------------------------------------------------------------
# 3. Expand genes that share the same genomic junction
# ------------------------------------------------------------

features["gene_list"] = (
    features["gene_id"]
    .astype(str)
    .str.split(";")
)

features_long = features.explode(
    "gene_list"
)

features_long["gene_list"] = (
    features_long["gene_list"]
    .str.strip()
)


# Rename gene column
features_long = features_long.rename(
    columns={
        "gene_list": "gene_id_master"
    }
)


# ------------------------------------------------------------
# 4. Load architecture-support table
# ------------------------------------------------------------

arch = pd.read_csv(
    ARCH_FILE,
    sep="\t"
)

print("Architecture rows:", len(arch))


# ------------------------------------------------------------
# 5. Check that gene-junction relationships match
# ------------------------------------------------------------

feature_pairs = set(
    zip(
        features_long["gene_id_master"],
        features_long["junction_id"]
    )
)

arch_pairs = set(
    zip(
        arch["gene_id"],
        arch["junction_id"]
    )
)


print(
    "Architecture pairs missing from features:",
    len(arch_pairs - feature_pairs)
)

print(
    "Feature pairs missing from architecture:",
    len(feature_pairs - arch_pairs)
)


# ------------------------------------------------------------
# 6. Keep only the columns needed for the master dataset
# ------------------------------------------------------------

feature_columns = [
    "gene_id_master",
    "gene_name",
    "junction_id",
    "intron_id",
    "chromosome",
    "strand",
    "intron_start",
    "intron_end",
    "donor_motif",
    "acceptor_motif",
    "splice_class",
    "transcripts",
    "gene_count",
    "exon_skipping",
    "alternative_donor",
    "alternative_acceptor",
    "structural_evidence"
]

features_clean = features_long[
    feature_columns
].copy()


# ------------------------------------------------------------
# 7. Rename gene ID to the final name
# ------------------------------------------------------------

features_clean = features_clean.rename(
    columns={
        "gene_id_master": "gene_id"
    }
)


# ------------------------------------------------------------
# 8. Merge architecture information
# ------------------------------------------------------------

master = features_clean.merge(
    arch,
    on=[
        "gene_id",
        "junction_id"
    ],
    how="left",
    validate="one_to_one"
)


# ------------------------------------------------------------
# 9. Check for missing architecture information
# ------------------------------------------------------------

missing_arch = master[
    master["architecture_fraction"].isna()
]

print(
    "Master rows:",
    len(master)
)

print(
    "Rows missing architecture information:",
    len(missing_arch)
)


# ------------------------------------------------------------
# 10. Remove duplicate rows if any
# ------------------------------------------------------------

master = master.drop_duplicates(
    subset=[
        "gene_id",
        "junction_id"
    ]
)


# ------------------------------------------------------------
# 11. Sort the dataset
# ------------------------------------------------------------

master = master.sort_values(
    [
        "chromosome",
        "intron_start",
        "intron_end",
        "gene_id"
    ]
)


# ------------------------------------------------------------
# 12. Save the master dataset
# ------------------------------------------------------------

master.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# 13. Final summary
# ------------------------------------------------------------

print("\n========================================")
print("MASTER DATASET CREATED")
print("========================================")

print("Output:", OUTPUT_FILE)
print("Rows:", len(master))
print("Unique genes:", master["gene_id"].nunique())
print("Unique junctions:", master["junction_id"].nunique())

print(
    "Rows with architecture fraction:",
    master["architecture_fraction"].notna().sum()
)

print(
    "Rows missing architecture fraction:",
    master["architecture_fraction"].isna().sum()
)

print("\nColumns:")
for column in master.columns:
    print(" -", column)
