import pandas as pd
import pyBigWig
import numpy as np


# ============================================================
# SCRIPT 18
# Full phyloP conservation extraction
# ============================================================

MASTER_FILE = "data/processed/master_gene_junction_table.tsv"

PHYLOP_FILE = "data/reference/conservation/hg38.phyloP100way.bw"

OUTPUT_FILE = "data/processed/junction_phyloP_conservation.tsv"


WINDOW = 20


# ------------------------------------------------------------
# 1. Load master dataset
# ------------------------------------------------------------

df = pd.read_csv(
    MASTER_FILE,
    sep="\t"
)

print("Master rows:", len(df))


# ------------------------------------------------------------
# 2. Open phyloP BigWig
# ------------------------------------------------------------

bw = pyBigWig.open(PHYLOP_FILE)

if bw is None:
    raise RuntimeError("Could not open phyloP BigWig file")

print("phyloP file opened successfully.")


# ------------------------------------------------------------
# 3. Function to obtain genomic scores
# ------------------------------------------------------------

def get_scores(chromosome, position):

    start = position - WINDOW - 1
    end = position + WINDOW

    return np.array(
        bw.values(
            chromosome,
            start,
            end
        ),
        dtype=float
    )


# ------------------------------------------------------------
# 4. Store results
# ------------------------------------------------------------

records = []


# ------------------------------------------------------------
# 5. Process every gene-junction relationship
# ------------------------------------------------------------

for i, row in df.iterrows():

    if i % 10000 == 0:
        print(
            f"Processed {i:,} / {len(df):,}"
        )

    chromosome = row["chromosome"]

    strand = row["strand"]

    junction = row["junction_id"]


    # --------------------------------------------------------
    # Parse junction coordinates
    # --------------------------------------------------------

    coordinates = junction.split(":")[1]

    donor = int(
        coordinates.split("-")[0]
    )

    acceptor = int(
        coordinates.split("-")[1]
    )


    # --------------------------------------------------------
    # Get donor and acceptor windows
    # --------------------------------------------------------

    donor_scores = get_scores(
        chromosome,
        donor
    )

    acceptor_scores = get_scores(
        chromosome,
        acceptor
    )


    # --------------------------------------------------------
    # Convert missing values to NaN
    # --------------------------------------------------------

    donor_valid = donor_scores[
        ~np.isnan(donor_scores)
    ]

    acceptor_valid = acceptor_scores[
        ~np.isnan(acceptor_scores)
    ]


    # --------------------------------------------------------
    # Summary statistics
    # --------------------------------------------------------

    if len(donor_valid) > 0:

        donor_mean = np.mean(donor_valid)
        donor_median = np.median(donor_valid)
        donor_min = np.min(donor_valid)
        donor_max = np.max(donor_valid)

    else:

        donor_mean = np.nan
        donor_median = np.nan
        donor_min = np.nan
        donor_max = np.nan


    if len(acceptor_valid) > 0:

        acceptor_mean = np.mean(acceptor_valid)
        acceptor_median = np.median(acceptor_valid)
        acceptor_min = np.min(acceptor_valid)
        acceptor_max = np.max(acceptor_valid)

    else:

        acceptor_mean = np.nan
        acceptor_median = np.nan
        acceptor_min = np.nan
        acceptor_max = np.nan


    # --------------------------------------------------------
    # Store complete profiles
    # --------------------------------------------------------

    donor_profile = ",".join(
        str(x) for x in donor_scores
    )

    acceptor_profile = ",".join(
        str(x) for x in acceptor_scores
    )


    # --------------------------------------------------------
    # Save record
    # --------------------------------------------------------

    records.append(
        {
            "gene_id": row["gene_id"],
            "gene_name": row["gene_name"],
            "junction_id": junction,
            "intron_id": row["intron_id"],
            "chromosome": chromosome,
            "strand": strand,

            "donor": donor,
            "acceptor": acceptor,

            "donor_motif": row["donor_motif"],
            "acceptor_motif": row["acceptor_motif"],

            "splice_class": row["splice_class"],

            "architecture_fraction":
                row["architecture_fraction"],

            "architecture_support":
                row["architecture_support"],

            "total_architectures":
                row["total_architectures"],

            "structural_evidence":
                row["structural_evidence"],

            "donor_mean_phyloP":
                donor_mean,

            "donor_median_phyloP":
                donor_median,

            "donor_min_phyloP":
                donor_min,

            "donor_max_phyloP":
                donor_max,

            "acceptor_mean_phyloP":
                acceptor_mean,

            "acceptor_median_phyloP":
                acceptor_median,

            "acceptor_min_phyloP":
                acceptor_min,

            "acceptor_max_phyloP":
                acceptor_max,

            "donor_valid_positions":
                len(donor_valid),

            "acceptor_valid_positions":
                len(acceptor_valid),

            "donor_profile":
                donor_profile,

            "acceptor_profile":
                acceptor_profile
        }
    )


# ------------------------------------------------------------
# 6. Close BigWig
# ------------------------------------------------------------

bw.close()


# ------------------------------------------------------------
# 7. Convert to DataFrame
# ------------------------------------------------------------

result = pd.DataFrame(records)


# ------------------------------------------------------------
# 8. Save
# ------------------------------------------------------------

result.to_csv(
    OUTPUT_FILE,
    sep="\t",
    index=False
)


# ------------------------------------------------------------
# 9. Final QC
# ------------------------------------------------------------

print("\n========================================")
print("FULL PHYLOP EXTRACTION COMPLETE")
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

print(
    "Rows with donor conservation:",
    result["donor_mean_phyloP"].notna().sum()
)

print(
    "Rows with acceptor conservation:",
    result["acceptor_mean_phyloP"].notna().sum()
)

print(
    "\nDonor mean phyloP summary:"
)

print(
    result["donor_mean_phyloP"].describe()
)

print(
    "\nAcceptor mean phyloP summary:"
)

print(
    result["acceptor_mean_phyloP"].describe()
)
