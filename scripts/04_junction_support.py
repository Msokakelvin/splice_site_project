import pandas as pd

junction_file = "data/processed/splice_junctions.tsv"
output_file = "data/processed/junction_support.tsv"

df = pd.read_csv(junction_file, sep="\t")

# Create a unique identifier for each genomic splice junction
df["junction_id"] = (
    df["chromosome"].astype(str)
    + ":"
    + df["donor"].astype(str)
    + "-"
    + df["acceptor"].astype(str)
    + ":"
    + df["strand"].astype(str)
)

# Number of transcripts belonging to each gene
gene_transcript_counts = (
    df[["gene_id", "transcript_id"]]
    .drop_duplicates()
    .groupby("gene_id")
    .size()
    .rename("gene_transcript_count")
)

# Number of distinct transcripts supporting each junction
junction_counts = (
    df[["gene_id", "junction_id", "transcript_id"]]
    .drop_duplicates()
    .groupby(["gene_id", "junction_id"])
    .size()
    .rename("supporting_transcripts")
)

# Combine the two measurements
support = junction_counts.reset_index()

support = support.merge(
    gene_transcript_counts.reset_index(),
    on="gene_id",
    how="left"
)

# Calculate proportion of transcripts supporting the junction
support["support_fraction"] = (
    support["supporting_transcripts"]
    / support["gene_transcript_count"]
)

# Classify the junction based on annotation support
support["support_class"] = "intermediate"

support.loc[
    support["supporting_transcripts"] == 1,
    "support_class"
] = "transcript_specific"

support.loc[
    support["supporting_transcripts"] == support["gene_transcript_count"],
    "support_class"
] = "shared_all_transcripts"

support.to_csv(output_file, sep="\t", index=False)

print("Unique junctions:", len(support))
print(
    "Genes represented:",
    support["gene_id"].nunique()
)
print(
    "Transcript-specific junctions:",
    (support["support_class"] == "transcript_specific").sum()
)
print(
    "Shared across all transcripts:",
    (support["support_class"] == "shared_all_transcripts").sum()
)
print("Output:", output_file)
