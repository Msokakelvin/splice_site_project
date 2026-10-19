import gzip
import pandas as pd

gtf_file = "data/reference/gencode.v50.annotation.gtf.gz"
output_file = "data/processed/transcripts.tsv"

records = []

with gzip.open(gtf_file, "rt") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if fields[2] != "transcript":
            continue

        attributes = fields[8]

        if 'gene_type "protein_coding"' not in attributes:
            continue

        gene_id = attributes.split('gene_id "')[1].split('"')[0]
        transcript_id = attributes.split('transcript_id "')[1].split('"')[0]
        gene_name = attributes.split('gene_name "')[1].split('"')[0]

        records.append([
            gene_id,
            gene_name,
            transcript_id,
            fields[0],
            int(fields[3]),
            int(fields[4]),
            fields[6]
        ])

df = pd.DataFrame(records, columns=[
    "gene_id",
    "gene_name",
    "transcript_id",
    "chromosome",
    "start",
    "end",
    "strand"
])

df.to_csv(output_file, sep="\t", index=False)

print("Transcript records:", len(df))
print("Unique genes:", df["gene_id"].nunique())
print("Output:", output_file)
