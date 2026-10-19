import pandas as pd

junction_file = "data/processed/splice_junctions.tsv"
output_file = "data/processed/intron_coordinates.tsv"

df = pd.read_csv(
    junction_file,
    sep="\t"
)

records = []

for transcript_id, group in df.groupby("transcript_id"):

    strand = group["strand"].iloc[0]

    if strand == "+":

        group = group.sort_values("upstream_exon_start")

    else:

        group = group.sort_values(
            "upstream_exon_start",
            ascending=False
        )

    for _, row in group.iterrows():

        if strand == "+":

            intron_start = int(row["upstream_exon_end"]) + 1
            intron_end = int(row["downstream_exon_start"]) - 1

        else:

            intron_start = int(row["downstream_exon_end"]) + 1
            intron_end = int(row["upstream_exon_start"]) - 1

        records.append([
            row["gene_id"],
            row["gene_name"],
            transcript_id,
            row["chromosome"],
            strand,
            intron_start,
            intron_end
        ])

result = pd.DataFrame(
    records,
    columns=[
        "gene_id",
        "gene_name",
        "transcript_id",
        "chromosome",
        "strand",
        "intron_start",
        "intron_end"
    ]
)

result.to_csv(
    output_file,
    sep="\t",
    index=False
)

print("Introns:", len(result))
print("Transcripts:", result["transcript_id"].nunique())
print("Output:", output_file)
