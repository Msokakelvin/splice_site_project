import pandas as pd
import numpy as np


# ============================================================
# SCRIPT 19
# Build biologically meaningful splice-site conservation
# features from the 41-position phyloP profiles
# ============================================================

INPUT_FILE = (
    "data/processed/junction_phyloP_conservation.tsv"
)

OUTPUT_FILE = (
    "data/processed/splice_conservation_features.tsv"
)


# ------------------------------------------------------------
# 1. Load data
# ------------------------------------------------------------

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)

print("Input rows:", len(df))


# ------------------------------------------------------------
# 2. Function to convert profile string to numbers
# ------------------------------------------------------------

def parse_profile(profile):

    if pd.isna(profile):
        return np.full(41, np.nan)

    values = [
        float(x)
        for x in str(profile).split(",")
    ]

    return np.array(values, dtype=float)


# ------------------------------------------------------------
# 3. Calculate features
# ------------------------------------------------------------

results = []


for i, row in df.iterrows():

    if i % 50000 == 0:
        print(
            f"Processing {i:,} / {len(df):,}"
        )

    strand = row["strand"]

    donor = parse_profile(
        row["donor_profile"]
    )

    acceptor = parse_profile(
        row["acceptor_profile"]
    )


    # --------------------------------------------------------
    # IMPORTANT:
    #
    # For negative-strand genes, genomic order is reversed.
    #
    # phyloP is not strand-specific, but the POSITIONAL
    # interpretation is strand-specific.
    #
    # Reverse the profiles so everything is represented
    # in transcript orientation.
    # --------------------------------------------------------

    if strand == "-":

        donor = donor[::-1]

        acceptor = acceptor[::-1]


    # --------------------------------------------------------
    # Position structure
    #
    # index 0  = -20
    # index 19 = -1
    # index 20 = 0
    # index 21 = +1
    # index 22 = +2
    # ...
    # index 40 = +20
    # --------------------------------------------------------


    # ========================================================
    # DONOR
    # ========================================================

    # upstream exon flank
    donor_exon = donor[0:21]

    # first 20 intronic bases
    donor_intron = donor[21:41]

    # GT splice motif
    donor_motif = donor[21:23]


    # ========================================================
    # ACCEPTOR
    # ========================================================

    # last 20 intronic bases
    acceptor_intron = acceptor[0:20]

    # AG splice motif
    acceptor_motif = acceptor[18:20]

    # downstream exon
    acceptor_exon = acceptor[20:41]


    # --------------------------------------------------------
    # Helper function
    # --------------------------------------------------------

    def mean_score(values):

        values = values[
            ~np.isnan(values)
        ]

        if len(values) == 0:
            return np.nan

        return np.mean(values)


    # --------------------------------------------------------
    # Calculate regional means
    # --------------------------------------------------------

    donor_exon_mean = mean_score(
        donor_exon
    )

    donor_intron_mean = mean_score(
        donor_intron
    )

    donor_motif_mean = mean_score(
        donor_motif
    )


    acceptor_intron_mean = mean_score(
        acceptor_intron
    )

    acceptor_motif_mean = mean_score(
        acceptor_motif
    )

    acceptor_exon_mean = mean_score(
        acceptor_exon
    )


    # --------------------------------------------------------
    # Overall splice-region means
    # --------------------------------------------------------

    donor_region_mean = mean_score(
        donor
    )

    acceptor_region_mean = mean_score(
        acceptor
    )


    # --------------------------------------------------------
    # Combined splice motif
    #
    # Donor GT + acceptor AG
    # --------------------------------------------------------

    splice_motif_values = np.concatenate(
        [
            donor_motif,
            acceptor_motif
        ]
    )

    splice_motif_mean = mean_score(
        splice_motif_values
    )


    # --------------------------------------------------------
    # Save result
    # --------------------------------------------------------

    results.append(
        {
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
            "junction_id": row["junction_id"],
            "intron_id": row["intron_id"],
            "chromosome": row["chromosome"],
            "strand": strand,

            "donor_motif":
                row["donor_motif"],

            "acceptor_motif":
                row["acceptor_motif"],

            "splice_class":
                row["splice_class"],

            "architecture_fraction":
                row["architecture_fraction"],

            "architecture_support":
                row["architecture_support"],

            "total_architectures":
                row["total_architectures"],

            "structural_evidence":
                row["structural_evidence"],

            # Donor
            "donor_exon_mean_phyloP":
                donor_exon_mean,

            "donor_motif_mean_phyloP":
                donor_motif_mean,

            "donor_intron_mean_phyloP":
                donor_intron_mean,

            "donor_region_mean_phyloP":
                donor_region_mean,

            # Acceptor
            "acceptor_intron_mean_phyloP":
                acceptor_intron_mean,

            "acceptor_motif_mean_phyloP":
                acceptor_motif_mean,

            "acceptor_exon_mean_phyloP":
                acceptor_exon_mean,

            "acceptor_region_mean_phyloP":
                acceptor_region_mean,

            # Combined
            "splice_motif_mean_phyloP":
                splice_motif_mean
        }
    )


# ------------------------------------------------------------
# 4. Convert to DataFrame
# ------------------------------------------------------------

result = pd.DataFrame(results)


# ------------------------------------------------------------
# 5. Save
# ------------------------------------------------------------

result.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# 6. QC
# ------------------------------------------------------------

print("\n========================================")
print("CONSERVATION FEATURES CREATED")
print("========================================")

print(
    "Output:",
    OUTPUT_FILE
)

print(
    "Rows:",
    len(result)
)

print(
    "Unique genes:",
    result["gene_id"].nunique()
)

print(
    "Unique junctions:",
    result["junction_id"].nunique()
)


print("\nMissing-value summary:")

columns_to_check = [
    "donor_exon_mean_phyloP",
    "donor_motif_mean_phyloP",
    "donor_intron_mean_phyloP",
    "acceptor_intron_mean_phyloP",
    "acceptor_motif_mean_phyloP",
    "acceptor_exon_mean_phyloP",
    "splice_motif_mean_phyloP"
]

for column in columns_to_check:

    missing = result[column].isna().sum()

    print(
        f"{column}: {missing:,} missing"
    )


print("\nSummary statistics:")

print(
    result[
        columns_to_check
    ].describe()
)
