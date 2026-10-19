import pandas as pd
import pyBigWig
import numpy as np

INPUT_FILE = "data/processed/master_gene_junction_table.tsv"
PHYLOP_FILE = "data/reference/conservation/hg38.phyloP100way.bw"
OUTPUT_FILE = "data/processed/correct_phyloP_profiles.tsv"

WINDOW = 20

df = pd.read_csv(INPUT_FILE, sep="\t")

bw = pyBigWig.open(PHYLOP_FILE)

print("Extracting correctly oriented phyloP profiles...")
print("Rows:", len(df))


def get_score(chrom, position):
    """
    Get phyloP score at a 1-based genomic position.
    """

    values = bw.values(
        chrom,
        position - 1,
        position
    )

    if len(values) == 0:
        return np.nan

    return values[0]


def donor_profile(row):
    """
    Return donor conservation in transcript orientation.

    Positions:
    -20 ... -1 | +1 +2 ... +20
                  G  T
    """

    chrom = row["chromosome"]
    strand = row["strand"]

    intron_start = int(row["intron_start"])
    intron_end = int(row["intron_end"])

    scores = []

    if strand == "+":

        # +1 = intron_start
        # +2 = intron_start + 1
        #
        # -1 = intron_start - 1
        # -2 = intron_start - 2
        #
        # Increasing genomic coordinates = transcript direction

        positions = range(
            intron_start - WINDOW,
            intron_start + WINDOW
        )

    else:

        # Negative strand:
        #
        # +1 = intron_end
        # +2 = intron_end - 1
        #
        # -1 = intron_end + 1
        # -2 = intron_end + 2
        #
        # Therefore genomic positions must be reversed.

        positions = range(
            intron_end + WINDOW,
            intron_end - WINDOW,
            -1
        )

    for pos in positions:
        scores.append(
            get_score(chrom, pos)
        )

    return scores


def acceptor_profile(row):
    """
    Return acceptor conservation in transcript orientation.

    Positions:
    -20 ... -2 -1 | +1 ... +20
                  A  G
    """

    chrom = row["chromosome"]
    strand = row["strand"]

    intron_start = int(row["intron_start"])
    intron_end = int(row["intron_end"])

    scores = []

    if strand == "+":

        # -2 = intron_end - 1
        # -1 = intron_end
        #
        # +1 = intron_end + 1

        positions = range(
            intron_end - WINDOW,
            intron_end + WINDOW
        )

    else:

        # Negative strand:
        #
        # -2 = intron_start + 1
        # -1 = intron_start
        #
        # +1 = intron_start - 1
        #
        # Transcript direction is decreasing genomic coordinate.

        positions = range(
            intron_start + WINDOW,
            intron_start - WINDOW,
            -1
        )

    for pos in positions:
        scores.append(
            get_score(chrom, pos)
        )

    return scores


records = []

for i, (_, row) in enumerate(df.iterrows()):

    donor = donor_profile(row)
    acceptor = acceptor_profile(row)

    records.append({
        "gene_id": row["gene_id"],
        "gene_name": row["gene_name"],
        "junction_id": row["junction_id"],
        "strand": row["strand"],
        "intron_start": row["intron_start"],
        "intron_end": row["intron_end"],
        "architecture_fraction": row["architecture_fraction"],
        "donor_profile": ",".join(
            "nan" if pd.isna(x) else f"{x:.6f}"
            for x in donor
        ),
        "acceptor_profile": ",".join(
            "nan" if pd.isna(x) else f"{x:.6f}"
            for x in acceptor
        )
    })

    if (i + 1) % 50000 == 0:
        print(
            f"Processed {i + 1:,} / {len(df):,}"
        )


result = pd.DataFrame(records)

result.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)

bw.close()

print("\n========================================")
print("CORRECT PHYLOP EXTRACTION COMPLETE")
print("========================================")

print("Output:", OUTPUT_FILE)
print("Rows:", len(result))
print("Profile length: 41 positions")

print("\nCoordinate system:")
print("Donor:    -20 ... -1 | +1 +2 ... +20")
print("Acceptor: -20 ... -2 -1 | +1 ... +20")
