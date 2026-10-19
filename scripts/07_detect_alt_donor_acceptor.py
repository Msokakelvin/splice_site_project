import pandas as pd

input_file = "data/processed/splice_junctions.tsv"
output_file = "data/processed/alternative_donor_acceptor_events.tsv"

# --------------------------------------------------
# 1. Load junctions
# --------------------------------------------------

df = pd.read_csv(input_file, sep="\t")


# --------------------------------------------------
# 2. Remove duplicate transcript-junction records
# --------------------------------------------------

df = df.drop_duplicates(
    subset=[
        "gene_id",
        "transcript_id",
        "donor",
        "acceptor"
    ]
)


# --------------------------------------------------
# 3. Alternative donor
#
# Same acceptor + different donors
# --------------------------------------------------

donor_groups = (
    df.groupby(
        ["gene_id", "chromosome", "strand", "acceptor"]
    )["donor"]
    .nunique()
    .reset_index(name="number_of_donors")
)

alternative_donor = donor_groups[
    donor_groups["number_of_donors"] > 1
].copy()

alternative_donor["event_type"] = "alternative_5prime_donor"


# --------------------------------------------------
# 4. Alternative acceptor
#
# Same donor + different acceptors
# --------------------------------------------------

acceptor_groups = (
    df.groupby(
        ["gene_id", "chromosome", "strand", "donor"]
    )["acceptor"]
    .nunique()
    .reset_index(name="number_of_acceptors")
)

alternative_acceptor = acceptor_groups[
    acceptor_groups["number_of_acceptors"] > 1
].copy()

alternative_acceptor["event_type"] = "alternative_3prime_acceptor"


# --------------------------------------------------
# 5. Save the two event types
# --------------------------------------------------

events = pd.concat(
    [
        alternative_donor,
        alternative_acceptor
    ],
    ignore_index=True
)

events.to_csv(
    output_file,
    sep="\t",
    index=False
)


# --------------------------------------------------
# 6. Summary
# --------------------------------------------------

print(
    "Alternative donor groups:",
    len(alternative_donor)
)

print(
    "Genes with alternative donor groups:",
    alternative_donor["gene_id"].nunique()
)

print(
    "Alternative acceptor groups:",
    len(alternative_acceptor)
)

print(
    "Genes with alternative acceptor groups:",
    alternative_acceptor["gene_id"].nunique()
)

print(
    "Total event groups:",
    len(events)
)

print(
    "Output:",
    output_file
)
