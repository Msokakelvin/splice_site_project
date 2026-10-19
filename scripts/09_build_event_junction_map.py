import pandas as pd

event_file = "data/processed/strict_exon_skipping_events.tsv"
junction_file = "data/processed/splice_junctions.tsv"

output_file = (
    "data/processed/exon_skipping_event_junctions.tsv"
)

print("Reading junction table...")

junctions = pd.read_csv(
    junction_file,
    sep="\t",
    usecols=[
        "gene_id",
        "transcript_id",
        "chromosome",
        "strand",
        "donor",
        "acceptor"
    ]
)

print("Junction rows:", len(junctions))


# --------------------------------------------------
# 1. Build transcript metadata lookup
# --------------------------------------------------

transcript_metadata = (
    junctions[
        [
            "transcript_id",
            "chromosome",
            "strand"
        ]
    ]
    .drop_duplicates()
    .set_index("transcript_id")
    .to_dict("index")
)

print(
    "Transcript metadata entries:",
    len(transcript_metadata)
)


# --------------------------------------------------
# 2. Read exon-skipping events
# --------------------------------------------------

print("Reading exon-skipping events...")

events = pd.read_csv(
    event_file,
    sep="\t"
)

events = events.drop_duplicates()

print("Events:", len(events))


# --------------------------------------------------
# 3. Build event-junction records
# --------------------------------------------------

records = []

for event_id, row in events.iterrows():

    gene_id = row["gene_id"]

    tx_included = row["transcript_with_exon"]
    tx_skipping = row["transcript_skipping_exon"]

    upstream_start = int(row["upstream_exon_start"])
    upstream_end = int(row["upstream_exon_end"])

    skipped_start = int(row["skipped_exon_start"])
    skipped_end = int(row["skipped_exon_end"])

    downstream_start = int(row["downstream_exon_start"])
    downstream_end = int(row["downstream_exon_end"])


    # ----------------------------------------------
    # Get transcript information
    # ----------------------------------------------

    metadata = transcript_metadata.get(
        tx_included
    )

    if metadata is None:
        continue

    chromosome = metadata["chromosome"]
    strand = metadata["strand"]


    # ----------------------------------------------
    # Determine junction coordinates
    # ----------------------------------------------

    if strand == "+":

        # Junction entering skipped exon
        j1 = (
            tx_included,
            upstream_end,
            skipped_start
        )

        # Junction leaving skipped exon
        j2 = (
            tx_included,
            skipped_end,
            downstream_start
        )

        # Junction created by exon skipping
        j3 = (
            tx_skipping,
            upstream_end,
            downstream_start
        )

    else:

        # Junction entering skipped exon
        j1 = (
            tx_included,
            upstream_start,
            skipped_end
        )

        # Junction leaving skipped exon
        j2 = (
            tx_included,
            skipped_start,
            downstream_end
        )

        # Junction created by exon skipping
        j3 = (
            tx_skipping,
            upstream_start,
            downstream_end
        )


    # ----------------------------------------------
    # Store the three event junctions
    # ----------------------------------------------

    event_junctions = [
        ("exon_entry", j1),
        ("exon_exit", j2),
        ("exon_skipping", j3)
    ]


    for role, junction in event_junctions:

        transcript_id = junction[0]
        donor = junction[1]
        acceptor = junction[2]

        records.append([
            gene_id,
            event_id,
            role,
            transcript_id,
            chromosome,
            strand,
            donor,
            acceptor
        ])


# --------------------------------------------------
# 4. Create output table
# --------------------------------------------------

result = pd.DataFrame(
    records,
    columns=[
        "gene_id",
        "event_id",
        "junction_role",
        "transcript_id",
        "chromosome",
        "strand",
        "donor",
        "acceptor"
    ]
)


# --------------------------------------------------
# 5. Save output
# --------------------------------------------------

result.to_csv(
    output_file,
    sep="\t",
    index=False
)


# --------------------------------------------------
# 6. Summary
# --------------------------------------------------

print()
print("====================================")
print("EVENT-JUNCTION MAP COMPLETE")
print("====================================")

print(
    "Event-junction rows:",
    len(result)
)

print(
    "Unique genes:",
    result["gene_id"].nunique()
)

print(
    "Unique junctions:",
    result[
        [
            "chromosome",
            "strand",
            "donor",
            "acceptor"
        ]
    ]
    .drop_duplicates()
    .shape[0]
)

print()
print("Junction roles:")

print(
    result["junction_role"].value_counts()
)

print()
print("Output:", output_file)
