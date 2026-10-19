import gzip
import pandas as pd

gtf_file = "data/reference/gencode.v50.annotation.gtf.gz"
output_file = "data/processed/exons.tsv"

records = []

with gzip.open(gtf_file, "rt") as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        # Keep only exon records
        if fields[2] != "exon":
            continue

        attributes = fields[8]

        # Keep only protein-coding genes
        if 'gene_type "protein_coding"' not in attributes:
            continue

        gene_id = attributes.split('gene_id "')[1].split('"')[0]
        transcript_id = attributes.split('transcript_id "')[1].split('"')[0]
        gene_name = attributes.split('gene_name "')[1].split('"')[0]

        exon_number = None

        if 'exon_number "' in attributes:
            exon_number = attributes.split('exon_number "')[1].split('"')[0]

        records.append([
            gene_id,
            gene_name,
            transcript_id,
            fields[0],
            int(fields[3]),
            int(fields[4]),
            fields[6],
            exon_number
        ])

df = pd.DataFrame(records, columns=[
    "gene_id",
    "gene_name",
    "transcript_id",
    "chromosome",
    "start",
    "end",
    "strand",
    "exon_number"
])

df.to_csv(output_file, sep="\t", index=False)

print("Exon records:", len(df))
print("Unique genes:", df["gene_id"].nunique())
print("Unique transcripts:", df["transcript_id"].nunique())
print("Output:", output_file)
