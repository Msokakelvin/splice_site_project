import pandas as pd

input_file = (
    "data/processed/exon_skipping_event_junctions.tsv"
)

output_file = (
    "data/processed/alternative_junction_summary.tsv"
)

print("Reading event-junction table...")

df = pd.read_csv(
    input_file,
    sep="\t"
)

print("Rows:", len(df))

# Create a unique identifier for each genomic splice junction
df["junction_id"] = (
    df["chromosome"].astype(str)
    + ":"
    + df["donor"].astype(int).astype(str)
    + "-"
    + df["acceptor"].astype(int).astype(str)
    + ":"
    + df["strand"].astype(str)
)

# Summarize each unique junction
summary = (
    df.groupby("junction_id")
    .agg(
        gene_id=("gene_id", lambda x: ";".join(sorted(set(x)))),
        chromosome=("chromosome", "first"),
        strand=("strand", "first"),
        donor=("donor", "first"),
        acceptor=("acceptor", "first"),
        events=("event_id", "nunique"),
        transcripts=("transcript_id", "nunique"),
        roles=("junction_role", lambda x: ";".join(sorted(set(x))))
    )
    .reset_index()
)

# Number of genes associated with each junction
summary["gene_count"] = (
    summary["gene_id"]
    .str.split(";")
    .str.len()
)

# Classify whether a junction participates in multiple roles
summary["role_class"] = "single_role"

summary.loc[
    summary["roles"].str.contains(";"),
    "role_class"
] = "multiple_roles"

# Save
summary.to_csv(
    output_file,
    sep="\t",
    index=False
)

print()
print("Unique junctions:", len(summary))
print(
    "Unique genes:",
    df["gene_id"].nunique()
)
print(
    "Junctions associated with multiple genes:",
    (summary["gene_count"] > 1).sum()
)
print(
    "Junctions with multiple roles:",
    (summary["role_class"] == "multiple_roles").sum()
)

print()
print("Role combinations:")
print(summary["roles"].value_counts())

print()
print("Output:", output_file)
