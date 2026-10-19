import pandas as pd

exon_file = "data/processed/exons.tsv"
output_file = "data/processed/exon_skipping_events.tsv"

df = pd.read_csv(exon_file, sep="\t")

events = []

# Process one gene at a time
for gene_id, gene_df in df.groupby("gene_id"):

    transcript_ids = gene_df["transcript_id"].unique()

    if len(transcript_ids) < 2:
        continue

    # Store exon coordinates for every transcript
    transcript_exons = {}

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

    # Compare pairs of transcripts
    for i in range(len(transcript_ids)):

        for j in range(i + 1, len(transcript_ids)):

            t1 = transcript_ids[i]
            t2 = transcript_ids[j]

            exons1 = transcript_exons[t1]
            exons2 = transcript_exons[t2]

            set1 = set(exons1)
            set2 = set(exons2)

            # Exons present in transcript 1 but absent in transcript 2
            unique_exons = set1 - set2

            for exon in unique_exons:

                index = exons1.index(exon)

                # The exon must be internal
                if index == 0 or index == len(exons1) - 1:
                    continue

                upstream_exon = exons1[index - 1]
                downstream_exon = exons1[index + 1]

                # Both neighboring exons must occur in transcript 2
                if (
                    upstream_exon in set2
                    and downstream_exon in set2
                ):

                    events.append([
                        gene_id,
                        t1,
                        t2,
                        exon[0],
                        exon[1],
                        upstream_exon[0],
                        upstream_exon[1],
                        downstream_exon[0],
                        downstream_exon[1],
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

print("Exon-skipping events:", len(events_df))
print(
    "Genes with exon-skipping events:",
    events_df["gene_id"].nunique()
)
print("Output:", output_file)
