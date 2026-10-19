import pandas as pd

junction_file = "data/processed/splice_junctions.tsv"
output_file = "data/processed/transcript_architecture.tsv"

df = pd.read_csv(junction_file, sep="\t")

# Create a unique junction identifier
df["junction_id"] = (
    df["chromosome"].astype(str)
    + ":"
    + df["donor"].astype(int).astype(str)
    + "-"
    + df["acceptor"].astype(int).astype(str)
    + ":"
    + df["strand"].astype(str)
)

architectures = []

for transcript_id, group in df.groupby("transcript_id"):

    strand = group["strand"].iloc[0]

    # Put junctions in biological transcription order
    if strand == "+":
        group = group.sort_values("donor")
    else:
        group = group.sort_values("donor", ascending=False)

    junction_chain = "|".join(group["junction_id"])

    architectures.append([
        group["gene_id"].iloc[0],
        group["gene_name"].iloc[0],
        transcript_id,
        strand,
        junction_chain,
        len(group)
    ])

architecture = pd.DataFrame(
    architectures,
    columns=[
        "gene_id",
        "gene_name",
        "transcript_id",
        "strand",
        "junction_chain",
        "junction_count"
    ]
)

# Number of distinct architectures within each gene
architecture_counts = (
    architecture
    .groupby("gene_id")["junction_chain"]
    .nunique()
    .rename("distinct_architectures")
    .reset_index()
)

architecture = architecture.merge(
    architecture_counts,
    on="gene_id",
    how="left"
)

architecture["architecture_class"] = "single_architecture"

architecture.loc[
    architecture["distinct_architectures"] > 1,
    "architecture_class"
] = "multiple_architectures"

architecture.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("Transcript architectures:", len(architecture))
print("Genes represented:", architecture["gene_id"].nunique())
print(
    "Genes with multiple architectures:",
    architecture.drop_duplicates("gene_id")["distinct_architectures"].gt(1).sum()
)
print(
    "Genes with single architecture:",
    architecture.drop_duplicates("gene_id")["distinct_architectures"].eq(1).sum()
)
print("Output:", output_file)
