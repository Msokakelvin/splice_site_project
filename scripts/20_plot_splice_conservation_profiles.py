import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


INPUT_FILE = "data/processed/junction_phyloP_conservation.tsv"

DONOR_OUTPUT = "results/figures/donor_conservation_profile.png"
ACCEPTOR_OUTPUT = "results/figures/acceptor_conservation_profile.png"


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)

print("Rows:", len(df))


# ------------------------------------------------------------
# Convert profile text into numeric arrays
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
# Store profiles
# ------------------------------------------------------------

donor_profiles = []
acceptor_profiles = []


# ------------------------------------------------------------
# Read every junction
# ------------------------------------------------------------

for i, row in df.iterrows():

    if i % 50000 == 0:
        print(
            f"Reading {i:,} / {len(df):,}"
        )

    donor = parse_profile(
        row["donor_profile"]
    )

    acceptor = parse_profile(
        row["acceptor_profile"]
    )

    # Reverse negative-strand profiles so that
    # all profiles are represented in transcript direction.

    if row["strand"] == "-":

        donor = donor[::-1]
        acceptor = acceptor[::-1]

    donor_profiles.append(donor)
    acceptor_profiles.append(acceptor)


# ------------------------------------------------------------
# Convert lists to matrices
# ------------------------------------------------------------

donor_matrix = np.array(
    donor_profiles
)

acceptor_matrix = np.array(
    acceptor_profiles
)


print(
    "Donor matrix:",
    donor_matrix.shape
)

print(
    "Acceptor matrix:",
    acceptor_matrix.shape
)


# ------------------------------------------------------------
# Calculate mean conservation at each position
# ------------------------------------------------------------

donor_mean = np.nanmean(
    donor_matrix,
    axis=0
)

acceptor_mean = np.nanmean(
    acceptor_matrix,
    axis=0
)


# ------------------------------------------------------------
# Calculate standard error
# ------------------------------------------------------------

donor_n = np.sum(
    ~np.isnan(donor_matrix),
    axis=0
)

acceptor_n = np.sum(
    ~np.isnan(acceptor_matrix),
    axis=0
)


donor_sd = np.nanstd(
    donor_matrix,
    axis=0
)

acceptor_sd = np.nanstd(
    acceptor_matrix,
    axis=0
)


donor_se = donor_sd / np.sqrt(
    donor_n
)

acceptor_se = acceptor_sd / np.sqrt(
    acceptor_n
)


# ------------------------------------------------------------
# Position labels
# ------------------------------------------------------------

positions = np.arange(
    -20,
    21
)


# ============================================================
# DONOR PLOT
# ============================================================

plt.figure(
    figsize=(10, 6)
)

plt.plot(
    positions,
    donor_mean,
    linewidth=2
)

plt.fill_between(
    positions,
    donor_mean - donor_se,
    donor_mean + donor_se,
    alpha=0.2
)

plt.axvline(
    x=0,
    linestyle="--"
)

plt.axhline(
    y=0,
    linestyle=":"
)

plt.xlabel(
    "Position relative to donor boundary"
)

plt.ylabel(
    "Mean phyloP conservation"
)

plt.title(
    "Evolutionary Conservation Around Splice Donor Sites"
)

plt.tight_layout()

plt.savefig(
    DONOR_OUTPUT,
    dpi=300
)

plt.close()


# ============================================================
# ACCEPTOR PLOT
# ============================================================

plt.figure(
    figsize=(10, 6)
)

plt.plot(
    positions,
    acceptor_mean,
    linewidth=2
)

plt.fill_between(
    positions,
    acceptor_mean - acceptor_se,
    acceptor_mean + acceptor_se,
    alpha=0.2
)

plt.axvline(
    x=0,
    linestyle="--"
)

plt.axhline(
    y=0,
    linestyle=":"
)

plt.xlabel(
    "Position relative to acceptor boundary"
)

plt.ylabel(
    "Mean phyloP conservation"
)

plt.title(
    "Evolutionary Conservation Around Splice Acceptor Sites"
)

plt.tight_layout()

plt.savefig(
    ACCEPTOR_OUTPUT,
    dpi=300
)

plt.close()


# ------------------------------------------------------------
# Final output
# ------------------------------------------------------------

print()
print("========================================")
print("CONSERVATION PROFILE ANALYSIS COMPLETE")
print("========================================")

print(
    "Donor figure:",
    DONOR_OUTPUT
)

print(
    "Acceptor figure:",
    ACCEPTOR_OUTPUT
)

print()
print("Donor profile:")

for position, value in zip(
    positions,
    donor_mean
):

    print(
        f"{position:>3}: {value:.4f}"
    )


print()
print("Acceptor profile:")

for position, value in zip(
    positions,
    acceptor_mean
):

    print(
        f"{position:>3}: {value:.4f}"
    )
