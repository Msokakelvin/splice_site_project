import pandas as pd

junction_file = "data/processed/splice_junctions.tsv"
transcript_file = "data/processed/transcripts.tsv"
output_file = "data/processed/junction_support.tsv"

# --------------------------------------------------
# 1. Load the two datasets
# --------------------------------------------------

junctions = pd.read_csv(junction_file, sep="\t")
transcripts = pd.read_csv(transcript_file, sep="\t")


# --------------------------------------------------
# 2. Create a unique ID for each genomic splice junction
# --------------------------------------------------

junctions["junction_id"] = (
    junctions["chromosome"].astype(str)
    + ":"
    + junctions["donor"].astype(str)
    + "-"
    + junctions["acceptor"].astype(str)
    + ":"
    + junctions["strand"].astype(str)
)


# --------------------------------------------------
# 3. Count ALL protein-coding transcripts per gene
# --------------------------------------------------

gene_transcript_counts = (
    transcripts[["gene_id", "transcript_id"]]
    .drop_duplicates()
    .groupby("gene_id")
    .size()
    .rename("gene_transcript_count")
)


# --------------------------------------------------
# 4. Count transcripts supporting each junction
# --------------------------------------------------

junction_counts = (
    junctions[
        ["gene_id", "junction_id", "transcript_id"]
    ]
    .drop_duplicates()
    .groupby(["gene_id", "junction_id"])
    .size()
    .rename("supporting_transcripts")
)


# --------------------------------------------------
# 5. Convert counts into a table
# --------------------------------------------------

support = junction_counts.reset_index()


# --------------------------------------------------
# 6. Add total transcript count for each gene
# --------------------------------------------------

support = support.merge(
    gene_transcript_counts.reset_index(),
    on="gene_id",
    how="left"
)


# --------------------------------------------------
# 7. Calculate junction support fraction
# --------------------------------------------------

support["support_fraction"] = (
    support["supporting_transcripts"]
    / support["gene_transcript_count"]
)


# --------------------------------------------------
# 8. Classify junctions
# --------------------------------------------------

support["support_class"] = "intermediate"

support.loc[
    support["supporting_transcripts"] == 1,
    "support_class"
] = "transcript_specific"

support.loc[
    support["supporting_transcripts"]
    == support["gene_transcript_count"],
    "support_class"
] = "shared_all_transcripts"


# --------------------------------------------------
# 9. Save the result
# --------------------------------------------------

support.to_csv(
    output_file,
    sep="\t",
    index=False
)


# --------------------------------------------------
# 10. Print summary
# --------------------------------------------------

print("Unique junctions:", len(support))

print(
    "Genes represented:",
    support["gene_id"].nunique()
)

print(
    "Transcript-specific junctions:",
    (
        support["support_class"]
        == "transcript_specific"
    ).sum()
)

print(
    "Shared across all transcripts:",
    (
        support["support_class"]
        == "shared_all_transcripts"
    ).sum()
)

print(
    "Output:",
    output_file
)
