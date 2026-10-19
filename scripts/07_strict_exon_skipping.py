import pandas as pd

exon_file = "data/processed/exons.tsv"
junction_file = "data/processed/splice_junctions.tsv"
output_file = "data/processed/strict_exon_skipping_events.tsv"

exons = pd.read_csv(exon_file, sep="\t")
junctions = pd.read_csv(junction_file, sep="\t")

# Remove duplicate transcript-junction records
junctions = junctions.drop_duplicates(
    subset=["gene_id", "transcript_id", "donor", "acceptor"]
)

events = []

for gene_id, gene_df in exons.groupby("gene_id"):

    transcript_ids = gene_df["transcript_id"].unique()

    if len(transcript_ids) < 2:
        continue

    transcript_exons = {}
    transcript_junctions = {}

    # Build exon and junction structures for each transcript
    for transcript_id, group in gene_df.groupby("transcript_id"):

        strand = group["strand"].iloc[0]

        if strand == "+":
            group = group.sort_values("start")
        else:
            group = group.sort_values("start", ascending=False)

        transcript_exons[transcript_id] = [
            (int(row.start), int(row.end))
            for _, row in group.iterrows()
        ]

    gene_junctions = junctions[junctions["gene_id"] == gene_id]

    for transcript_id, group in gene_junctions.groupby("transcript_id"):

        transcript_junctions[transcript_id] = set(
            (int(row.donor), int(row.acceptor))
            for _, row in group.iterrows()
        )

    # Compare transcript pairs
    transcript_list = list(transcript_exons.keys())

    for i in range(len(transcript_list)):

        for j in range(i + 1, len(transcript_list)):

            t1 = transcript_list[i]
            t2 = transcript_list[j]

            exons1 = transcript_exons[t1]
            exons2 = transcript_exons[t2]

            junc1 = transcript_junctions.get(t1, set())
            junc2 = transcript_junctions.get(t2, set())

            set2 = set(exons2)

            # Search for an internal exon present in transcript 1
            # but absent from transcript 2
            for k in range(1, len(exons1) - 1):

                skipped = exons1[k]
                upstream = exons1[k - 1]
                downstream = exons1[k + 1]

                if skipped in set2:
                    continue

                if upstream not in set2 or downstream not in set2:
                    continue

                strand = gene_df["strand"].iloc[0]

                if strand == "+":

                    junction_up_to_skipped = (
                        upstream[1],
                        skipped[0]
                    )

                    junction_skipped_to_downstream = (
                        skipped[1],
                        downstream[0]
                    )

                    direct_junction = (
                        upstream[1],
                        downstream[0]
                    )

                else:

                    junction_up_to_skipped = (
                        upstream[0],
                        skipped[1]
                    )

                    junction_skipped_to_downstream = (
                        skipped[0],
                        downstream[1]
                    )

                    direct_junction = (
                        upstream[0],
                        downstream[1]
                    )

                # Transcript 1 must contain the two junctions
                # surrounding the exon
                included_structure = (
                    junction_up_to_skipped in junc1
                    and junction_skipped_to_downstream in junc1
                )

                # Transcript 2 must contain the direct junction
                skipping_structure = direct_junction in junc2

                if included_structure and skipping_structure:

                    events.append([
                        gene_id,
                        t1,
                        t2,
                        skipped[0],
                        skipped[1],
                        upstream[0],
                        upstream[1],
                        downstream[0],
                        downstream[1],
                        "exon_skipping"
                    ])

events_df = pd.DataFrame(
    events,
    columns=[
        "gene_id",
        "transcript_with_exon",
        "transcript_skipping_exon",
        "skipped_exon_start",
        "skipped_exon_end",
        "upstream_exon_start",
        "upstream_exon_end",
        "downstream_exon_start",
        "downstream_exon_end",
        "event_type"
    ]
)

events_df = events_df.drop_duplicates()

events_df.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("Strict transcript-pair events:", len(events_df))

if len(events_df) > 0:

    unique_events = events_df.drop_duplicates(
        subset=[
            "gene_id",
            "skipped_exon_start",
            "skipped_exon_end",
            "upstream_exon_start",
            "upstream_exon_end",
            "downstream_exon_start",
            "downstream_exon_end"
        ]
    )

    print("Strict unique exon-skipping events:", len(unique_events))
    print(
        "Genes with strict exon-skipping events:",
        unique_events["gene_id"].nunique()
    )

print("Output:", output_file)
