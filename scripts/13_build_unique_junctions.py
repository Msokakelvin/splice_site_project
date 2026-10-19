import pandas as pd

input_file = "data/processed/splice_site_sequences.tsv"
output_file = "data/processed/unique_splice_sites.tsv"

print("Reading splice-site dataset...")

df = pd.read_csv(
    input_file,
    sep="\t"
)

print("Transcript-specific introns:", len(df))


# Create genomic intron ID
df["intron_id"] = (
    df["chromosome"].astype(str)
    + ":"
    + df["intron_start"].astype(str)
    + "-"
    + df["intron_end"].astype(str)
    + ":"
    + df["strand"].astype(str)
)


# Collapse transcript redundancy
unique = (
    df.groupby("intron_id")
    .agg(
        gene_id=("gene_id", lambda x: ";".join(sorted(set(x)))),
        gene_name=("gene_name", lambda x: ";".join(sorted(set(x)))),
        chromosome=("chromosome", "first"),
        strand=("strand", "first"),
        intron_start=("intron_start", "first"),
        intron_end=("intron_end", "first"),
        donor_motif=("donor_motif", "first"),
        acceptor_motif=("acceptor_motif", "first"),
        splice_class=("splice_class", "first"),
        transcripts=("transcript_id", "nunique"),
    )
    .reset_index()
)


unique["gene_count"] = (
    unique["gene_id"]
    .str.split(";")
    .str.len()
)


unique.to_csv(
    output_file,
    sep="\t",
    index=False
)


print()
print("Unique splice sites:", len(unique))

print(
    "Genes represented:",
    df["gene_id"].nunique()
)

print(
    "Sites shared by multiple transcripts:",
    (unique["transcripts"] > 1).sum()
)

print(
    "Sites associated with multiple genes:",
    (unique["gene_count"] > 1).sum()
)

print()
print("Splice classes:")

print(
    unique["splice_class"].value_counts()
)

print()
print("Output:", output_file)
