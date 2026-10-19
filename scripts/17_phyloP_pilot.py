import pandas as pd
import pyBigWig


# ============================================================
# SCRIPT 17
# phyloP conservation pilot
# ============================================================

MASTER_FILE = "data/processed/master_gene_junction_table.tsv"

PHYLOP_FILE = "data/reference/conservation/hg38.phyloP100way.bw"

OUTPUT_FILE = "data/processed/phyloP_pilot.tsv"


# ------------------------------------------------------------
# 1. Load master dataset
# ------------------------------------------------------------

df = pd.read_csv(
    MASTER_FILE,
    sep="\t"
)

print("Master rows:", len(df))


# ------------------------------------------------------------
# 2. Select a small pilot set
# ------------------------------------------------------------

# We deliberately select different structural classes.

high = df[
    (df["architecture_fraction"] >= 0.90) &
    (df["splice_class"] == "canonical_GT_AG")
].head(3)

low = df[
    (df["architecture_fraction"] <= 0.10) &
    (df["splice_class"] == "canonical_GT_AG")
].head(3)

negative = df[
    (df["strand"] == "-") &
    (df["splice_class"] == "canonical_GT_AG")
].head(2)

noncanonical = df[
    df["splice_class"].isin(
        ["noncanonical_GC_AG", "noncanonical_AT_AC"]
    )
].head(2)


pilot = pd.concat(
    [high, low, negative, noncanonical]
).drop_duplicates(
    subset=["gene_id", "junction_id"]
)

print("Pilot junctions:", len(pilot))

print(
    pilot[
        [
            "gene_name",
            "junction_id",
            "strand",
            "donor_motif",
            "acceptor_motif",
            "architecture_fraction",
            "splice_class"
        ]
    ].to_string(index=False)
)


# ------------------------------------------------------------
# 3. Open phyloP file
# ------------------------------------------------------------

bw = pyBigWig.open(PHYLOP_FILE)

if bw is None:
    raise RuntimeError("Could not open phyloP BigWig file")

print("\nphyloP BigWig opened successfully.")


# ------------------------------------------------------------
# 4. Function to obtain conservation scores
# ------------------------------------------------------------

WINDOW = 20


def get_scores(chromosome, position):

    # BigWig uses 0-based, half-open coordinates.
    #
    # Biological position is 1-based.
    #
    # Therefore:
    # position - 20 through position + 20
    # becomes:
    #
    # start = position - 21
    # end   = position + 20

    start = position - WINDOW - 1
    end = position + WINDOW

    scores = bw.values(
        chromosome,
        start,
        end
    )

    return scores


# ------------------------------------------------------------
# 5. Extract donor and acceptor conservation
# ------------------------------------------------------------

records = []


for _, row in pilot.iterrows():

    chromosome = row["chromosome"]
    strand = row["strand"]

    junction = row["junction_id"]

    # Parse junction:
    #
    # chr1:65433-65520:+

    coordinates = junction.split(":")[1]

    donor = int(
        coordinates.split("-")[0]
    )

    acceptor = int(
        coordinates.split("-")[1]
    )


    # --------------------------------------------------------
    # Donor
    # --------------------------------------------------------

    donor_scores = get_scores(
        chromosome,
        donor
    )


    # --------------------------------------------------------
    # Acceptor
    # --------------------------------------------------------

    acceptor_scores = get_scores(
        chromosome,
        acceptor
    )


    records.append(
        {
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
            "junction_id": junction,
            "strand": strand,
            "donor": donor,
            "acceptor": acceptor,
            "donor_motif": row["donor_motif"],
            "acceptor_motif": row["acceptor_motif"],
            "architecture_fraction": row[
                "architecture_fraction"
            ],
            "splice_class": row["splice_class"],

            "donor_scores": ",".join(
                str(x) for x in donor_scores
            ),

            "acceptor_scores": ",".join(
                str(x) for x in acceptor_scores
            ),

            "donor_mean_phyloP": (
                sum(
                    x for x in donor_scores
                    if x is not None
                )
                /
                sum(
                    x is not None
                    for x in donor_scores
                )
            ),

            "acceptor_mean_phyloP": (
                sum(
                    x for x in acceptor_scores
                    if x is not None
                )
                /
                sum(
                    x is not None
                    for x in acceptor_scores
                )
            )
        }
    )


# ------------------------------------------------------------
# 6. Close BigWig
# ------------------------------------------------------------

bw.close()


# ------------------------------------------------------------
# 7. Save results
# ------------------------------------------------------------

result = pd.DataFrame(records)

result.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# 8. Final summary
# ------------------------------------------------------------

print("\n========================================")
print("PHYLOP PILOT COMPLETE")
print("========================================")

print("Output:", OUTPUT_FILE)
print("Rows:", len(result))

print(
    "\nConservation summary:"
)

print(
    result[
        [
            "gene_name",
            "junction_id",
            "strand",
            "architecture_fraction",
            "donor_mean_phyloP",
            "acceptor_mean_phyloP"
        ]
    ].to_string(index=False)
)
