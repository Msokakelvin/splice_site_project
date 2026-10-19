import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path


# ============================================================
# FILES
# ============================================================

INPUT_FILE = (
    "data/processed/correct_phyloP_profiles.tsv"
)

OUTPUT_TABLE = (
    "results/tables/splice_conservation_profiles.tsv"
)

OUTPUT_DONOR = (
    "results/figures/donor_conservation_profile.png"
)

OUTPUT_ACCEPTOR = (
    "results/figures/acceptor_conservation_profile.png"
)

OUTPUT_COMBINED = (
    "results/figures/donor_acceptor_conservation_profiles.png"
)


# ============================================================
# SETTINGS
# ============================================================

POSITIONS = list(range(-20, 0)) + list(range(1, 21))

WINDOW_LABEL = (
    "−20 ... −1 | +1 ... +20"
)


# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

Path("results/tables").mkdir(
    parents=True,
    exist_ok=True
)

Path("results/figures").mkdir(
    parents=True,
    exist_ok=True
)


# ============================================================
# LOAD DATA
# ============================================================

print()
print("=" * 70)
print("SCRIPT 25: SPLICE-SITE CONSERVATION PROFILES")
print("=" * 70)

print()
print("Loading:")
print(INPUT_FILE)

df = pd.read_csv(
    INPUT_FILE,
    sep="\t"
)

print(
    f"Rows loaded: {len(df):,}"
)


# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns = [
    "donor_profile",
    "acceptor_profile"
]

missing = [
    col
    for col in required_columns
    if col not in df.columns
]

if missing:

    raise ValueError(
        "Missing required columns: "
        + ", ".join(missing)
    )


# ============================================================
# PROFILE PARSER
# ============================================================

def parse_profile(value):

    if pd.isna(value):
        return None

    try:

        values = [
            float(x)
            for x in str(value).split(",")
        ]

        return values

    except ValueError:

        return None


# ============================================================
# PARSE DONOR AND ACCEPTOR PROFILES
# ============================================================

print()
print("Parsing profiles...")

df["donor_values"] = df["donor_profile"].apply(
    parse_profile
)

df["acceptor_values"] = df["acceptor_profile"].apply(
    parse_profile
)


# ============================================================
# CHECK PROFILE LENGTHS
# ============================================================

donor_lengths = (
    df["donor_values"]
    .dropna()
    .apply(len)
)

acceptor_lengths = (
    df["acceptor_values"]
    .dropna()
    .apply(len)
)

print()
print("Profile length check")
print("-" * 70)

print(
    "Donor lengths:"
)

print(
    donor_lengths.value_counts()
    .sort_index()
    .to_string()
)

print()

print(
    "Acceptor lengths:"
)

print(
    acceptor_lengths.value_counts()
    .sort_index()
    .to_string()
)


# ============================================================
# REQUIRE 40 POSITIONS
# ============================================================

if not donor_lengths.empty:

    if not (donor_lengths == 40).all():

        raise ValueError(
            "Some donor profiles do not contain exactly 40 values."
        )


if not acceptor_lengths.empty:

    if not (acceptor_lengths == 40).all():

        raise ValueError(
            "Some acceptor profiles do not contain exactly 40 values."
        )


# ============================================================
# PROFILE SUMMARY FUNCTION
# ============================================================

def summarize_profiles(
    profile_column,
    profile_type
):

    records = []

    for index, position in enumerate(POSITIONS):

        values = []

        for profile in df[profile_column]:

            if profile is None:
                continue

            if len(profile) != 40:
                continue

            value = profile[index]

            if np.isfinite(value):

                values.append(value)


        values = np.array(values)

        n = len(values)


        if n > 0:

            mean = np.mean(values)

            sd = np.std(
                values,
                ddof=1
            ) if n > 1 else np.nan

            se = (
                sd / np.sqrt(n)
                if n > 1
                else np.nan
            )

            median = np.median(values)

        else:

            mean = np.nan
            sd = np.nan
            se = np.nan
            median = np.nan


        records.append({

            "profile": profile_type,

            "position": position,

            "mean_phyloP": mean,

            "median_phyloP": median,

            "sd_phyloP": sd,

            "se_phyloP": se,

            "n": n

        })


    return pd.DataFrame(records)


# ============================================================
# BUILD SUMMARIES
# ============================================================

print()
print("Calculating donor profile...")

donor_summary = summarize_profiles(
    "donor_values",
    "donor"
)

print(
    "Calculating acceptor profile..."
)

acceptor_summary = summarize_profiles(
    "acceptor_values",
    "acceptor"
)


# ============================================================
# COMBINE
# ============================================================

summary = pd.concat(
    [
        donor_summary,
        acceptor_summary
    ],
    ignore_index=True
)


# ============================================================
# SAVE TABLE
# ============================================================

summary.to_csv(
    OUTPUT_TABLE,
    sep="\t",
    index=False
)

print()
print(
    f"Saved table: {OUTPUT_TABLE}"
)


# ============================================================
# PLOT FUNCTION
# ============================================================

def make_profile_plot(
    data,
    profile_name,
    output_file
):

    x = data["position"].values

    y = data["mean_phyloP"].values

    se = data["se_phyloP"].values


    plt.figure(
        figsize=(11, 6)
    )


    plt.plot(
        x,
        y,
        linewidth=2,
        label=profile_name
    )


    valid_se = np.isfinite(se)

    if valid_se.any():

        plt.fill_between(
            x[valid_se],
            (
                y[valid_se]
                -
                se[valid_se]
            ),
            (
                y[valid_se]
                +
                se[valid_se]
            ),
            alpha=0.25
        )


    # Splice boundary

    plt.axvline(
        x=0,
        linestyle="--",
        linewidth=1
    )


    # Motif region

    if profile_name.lower() == "donor":

        plt.axvspan(
            1,
            2,
            alpha=0.10
        )

    else:

        plt.axvspan(
            -2,
            -1,
            alpha=0.10
        )


    plt.xlabel(
        "Position relative to splice junction"
    )

    plt.ylabel(
        "Mean phyloP conservation score"
    )

    plt.title(
        f"{profile_name.capitalize()} splice-site conservation"
    )


    plt.xticks(
        list(range(-20, 21, 5))
    )


    plt.grid(
        alpha=0.25
    )


    plt.tight_layout()

    plt.savefig(
        output_file,
        dpi=300
    )

    plt.close()


    print(
        f"Saved figure: {output_file}"
    )


# ============================================================
# DONOR FIGURE
# ============================================================

make_profile_plot(
    donor_summary,
    "donor",
    OUTPUT_DONOR
)


# ============================================================
# ACCEPTOR FIGURE
# ============================================================

make_profile_plot(
    acceptor_summary,
    "acceptor",
    OUTPUT_ACCEPTOR
)


# ============================================================
# COMBINED FIGURE
# ============================================================

plt.figure(
    figsize=(12, 7)
)

x = donor_summary["position"].values


# Donor

donor_y = (
    donor_summary["mean_phyloP"].values
)

donor_se = (
    donor_summary["se_phyloP"].values
)


plt.plot(
    x,
    donor_y,
    linewidth=2,
    label="Donor"
)


valid = np.isfinite(donor_se)

if valid.any():

    plt.fill_between(
        x[valid],
        donor_y[valid] - donor_se[valid],
        donor_y[valid] + donor_se[valid],
        alpha=0.20
    )


# Acceptor

acceptor_y = (
    acceptor_summary["mean_phyloP"].values
)

acceptor_se = (
    acceptor_summary["se_phyloP"].values
)


plt.plot(
    x,
    acceptor_y,
    linewidth=2,
    label="Acceptor"
)


valid = np.isfinite(acceptor_se)

if valid.any():

    plt.fill_between(
        x[valid],
        acceptor_y[valid] - acceptor_se[valid],
        acceptor_y[valid] + acceptor_se[valid],
        alpha=0.20
    )


# Boundary

plt.axvline(
    x=0,
    linestyle="--",
    linewidth=1
)


plt.xlabel(
    "Position relative to splice junction"
)

plt.ylabel(
    "Mean phyloP conservation score"
)

plt.title(
    "Genome-wide conservation around splice junctions"
)

plt.xticks(
    list(range(-20, 21, 5))
)

plt.legend()

plt.grid(
    alpha=0.25
)

plt.tight_layout()


plt.savefig(
    OUTPUT_COMBINED,
    dpi=300
)

plt.close()


print(
    f"Saved figure: {OUTPUT_COMBINED}"
)


# ============================================================
# FINAL SUMMARY
# ============================================================

print()
print("=" * 70)
print("SCRIPT 25 COMPLETE")
print("=" * 70)

print()
print(
    f"Junction records analyzed: {len(df):,}"
)

print(
    f"Profile positions: {len(POSITIONS)}"
)

print(
    "Positions:"
)

print(
    POSITIONS
)

print()
print(
    "Output table:"
)

print(
    OUTPUT_TABLE
)

print()
print(
    "Figures:"
)

print(
    OUTPUT_DONOR
)

print(
    OUTPUT_ACCEPTOR
)

print(
    OUTPUT_COMBINED
)

print()
print("Done.")
