import pandas as pd

exon_file = "data/processed/exons.tsv"
output_file = "data/processed/splice_junctions.tsv"

df = pd.read_csv(exon_file, sep="\t")

junctions = []

for transcript_id, group in df.groupby("transcript_id"):

    # A transcript needs at least two exons
    if len(group) < 2:
        continue

    strand = group["strand"].iloc[0]

    # Put exons in transcript order
    if strand == "+":
        group = group.sort_values("start")
    else:
        group = group.sort_values("start", ascending=False)

    exons = group.to_dict("records")

    for i in range(len(exons) - 1):

        exon1 = exons[i]
        exon2 = exons[i + 1]

        if strand == "+":

            donor = exon1["end"]
            acceptor = exon2["start"]

        else:

            donor = exon1["start"]
            acceptor = exon2["end"]

        junctions.append([
            exon1["gene_id"],
            exon1["gene_name"],
            transcript_id,
            exon1["chromosome"],
            strand,
            exon1["start"],
            exon1["end"],
            exon2["start"],
            exon2["end"],
            donor,
            acceptor
        ])

junction_df = pd.DataFrame(junctions, columns=[
    "gene_id",
    "gene_name",
    "transcript_id",
    "chromosome",
    "strand",
    "upstream_exon_start",
    "upstream_exon_end",
    "downstream_exon_start",
    "downstream_exon_end",
    "donor",
    "acceptor"
])

junction_df.to_csv(output_file, sep="\t", index=False)

print("Splice junctions:", len(junction_df))
print("Unique genes:", junction_df["gene_id"].nunique())
print("Unique transcripts:", junction_df["transcript_id"].nunique())
print("Output:", output_file)
