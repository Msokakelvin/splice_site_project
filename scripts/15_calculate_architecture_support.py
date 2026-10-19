import pandas as pd


architecture_file = "data/processed/transcript_architecture.tsv"
output_file = "data/processed/splice_site_architecture_support.tsv"


print("Reading transcript architectures...")

architecture = pd.read_csv(
    architecture_file,
    sep="\t"
)

print("Transcript records:", len(architecture))
print(
    "Unique transcripts:",
    architecture["transcript_id"].nunique()
)
print(
    "Unique junction chains:",
    architecture["junction_chain"].nunique()
)


# ==========================================================
# 1. IDENTIFY DISTINCT ARCHITECTURES
# ==========================================================

architectures = (
    architecture[
        [
            "gene_id",
            "gene_name",
            "junction_chain"
        ]
    ]
    .drop_duplicates()
    .copy()
)

print(
    "Unique gene-level architectures:",
    len(architectures)
)


# ==========================================================
# 2. CREATE AN ARCHITECTURE ID
# ==========================================================

architectures["architecture_id"] = (
    architectures["gene_id"]
    + "|"
    + architectures["junction_chain"]
)


# ==========================================================
# 3. COUNT TOTAL ARCHITECTURES PER GENE
# ==========================================================

architecture_counts = (
    architectures
    .groupby("gene_id")
    .size()
    .rename("total_architectures")
    .reset_index()
)


# ==========================================================
# 4. SPLIT JUNCTION CHAINS
# ==========================================================

architectures["junction_id"] = (
    architectures["junction_chain"]
    .str.split("|")
)


# ==========================================================
# 5. EXPAND ARCHITECTURES INTO JUNCTIONS
# ==========================================================

expanded = architectures[
    [
        "gene_id",
        "gene_name",
        "architecture_id",
        "junction_id"
    ]
].explode(
    "junction_id"
)


expanded = expanded[
    expanded["junction_id"].notna()
    & (expanded["junction_id"] != "")
].copy()


# ==========================================================
# 6. REMOVE DUPLICATE JUNCTIONS WITHIN AN ARCHITECTURE
# ==========================================================

expanded = expanded.drop_duplicates(
    subset=[
        "gene_id",
        "architecture_id",
        "junction_id"
    ]
)


# ==========================================================
# 7. COUNT DISTINCT ARCHITECTURES
#    CONTAINING EACH JUNCTION
# ==========================================================

support = (
    expanded
    .groupby(
        [
            "gene_id",
            "junction_id"
        ]
    )["architecture_id"]
    .nunique()
    .rename(
        "architecture_support"
    )
    .reset_index()
)


# ==========================================================
# 8. ADD TOTAL ARCHITECTURES PER GENE
# ==========================================================

support = support.merge(
    architecture_counts,
    on="gene_id",
    how="left"
)


# ==========================================================
# 9. CALCULATE ARCHITECTURE FRACTION
# ==========================================================

support["architecture_fraction"] = (
    support["architecture_support"]
    /
    support["total_architectures"]
)


# ==========================================================
# 10. SAVE
# ==========================================================

support.to_csv(
    output_file,
    sep="\t",
    index=False
)


# ==========================================================
# 11. SUMMARY
# ==========================================================

print()
print("====================================")
print("ARCHITECTURE SUPPORT COMPLETE")
print("====================================")

print(
    "Architecture-junction records:",
    len(support)
)

print(
    "Genes represented:",
    support["gene_id"].nunique()
)

print()
print("Architecture support summary:")

print(
    support[
        "architecture_support"
    ].describe()
)

print()
print("Architecture fraction summary:")

print(
    support[
        "architecture_fraction"
    ].describe()
)

print()
print("Maximum architecture support:")

print(
    support[
        "architecture_support"
    ].max()
)

print()
print("Output:", output_file)
