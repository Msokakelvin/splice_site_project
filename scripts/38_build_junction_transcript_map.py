import pandas as pd

INPUT = "data/processed/splice_junctions.tsv"
OUTPUT = "data/processed/junction_transcript_map.tsv"

print("Reading transcript-specific junctions...")

usecols = [
    "gene_id",
    "transcript_id",
    "chromosome",
    "strand",
    "donor",
    "acceptor"
]

chunksize = 100000

header_written = False
total = 0

for chunk in pd.read_csv(
    INPUT,
    sep="\t",
    usecols=usecols,
    chunksize=chunksize
):

    # Create the same junction ID used by the master table
    chunk["junction_id"] = (
        chunk["chromosome"].astype(str)
        + ":"
        + chunk[["donor", "acceptor"]].min(axis=1).astype(str)
        + "-"
        + chunk[["donor", "acceptor"]].max(axis=1).astype(str)
        + ":"
        + chunk["strand"].astype(str)
    )

    result = chunk[
        [
            "gene_id",
            "transcript_id",
            "junction_id"
        ]
    ].drop_duplicates()

    result.to_csv(
        OUTPUT,
        sep="\t",
        index=False,
        mode="w" if not header_written else "a",
        header=not header_written
    )

    header_written = True

    total += len(chunk)

    print(
        "Processed:",
        total
    )

print()
print("==============================================")
print("JUNCTION-TRANSCRIPT MAP COMPLETE")
print("==============================================")

print()
print("Output:")
print(OUTPUT)
