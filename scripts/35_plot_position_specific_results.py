import pandas as pd
import matplotlib.pyplot as plt

INPUT = "results/tables/position_specific_architecture_conservation.tsv"

DONOR_OUT = "results/figures/donor_position_architecture_effect.png"
ACCEPTOR_OUT = "results/figures/acceptor_position_architecture_effect.png"
COMBINED_OUT = "results/figures/position_specific_architecture_effects.png"


# ---------------------------------------------------------
# Load results
# ---------------------------------------------------------

df = pd.read_csv(INPUT, sep="\t")

print("Rows:", len(df))


# ---------------------------------------------------------
# Separate donor and acceptor
# ---------------------------------------------------------

donor = df[df["site"] == "donor"].copy()
acceptor = df[df["site"] == "acceptor"].copy()


# ---------------------------------------------------------
# Plot function
# ---------------------------------------------------------

def make_plot(data, title, output):

    data = data.sort_values("position")

    x = data["position"]

    y = data["coefficient"]

    lower = data["coefficient"] - data["CI_low"]

    upper = data["CI_high"] - data["coefficient"]

    plt.figure(figsize=(12, 6))

    plt.errorbar(
        x,
        y,
        yerr=[lower, upper],
        fmt="o-",
        capsize=3
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
        "Position relative to splice junction"
    )

    plt.ylabel(
        "Architecture-fraction coefficient"
    )

    plt.title(title)

    plt.xticks(
        range(-20, 21, 2)
    )

    plt.tight_layout()

    plt.savefig(
        output,
        dpi=300
    )

    plt.close()

    print("Saved:", output)


# ---------------------------------------------------------
# Donor
# ---------------------------------------------------------

make_plot(
    donor,
    "Donor: architecture-conservation association by position",
    DONOR_OUT
)


# ---------------------------------------------------------
# Acceptor
# ---------------------------------------------------------

make_plot(
    acceptor,
    "Acceptor: architecture-conservation association by position",
    ACCEPTOR_OUT
)


# ---------------------------------------------------------
# Combined figure
# ---------------------------------------------------------

plt.figure(figsize=(12, 6))

plt.errorbar(
    donor["position"],
    donor["coefficient"],
    yerr=[
        donor["coefficient"] - donor["CI_low"],
        donor["CI_high"] - donor["coefficient"]
    ],
    fmt="o-",
    capsize=3,
    label="Donor"
)

plt.errorbar(
    acceptor["position"],
    acceptor["coefficient"],
    yerr=[
        acceptor["coefficient"] - acceptor["CI_low"],
        acceptor["CI_high"] - acceptor["coefficient"]
    ],
    fmt="o-",
    capsize=3,
    label="Acceptor"
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
    "Position relative to splice junction"
)

plt.ylabel(
    "Architecture-fraction coefficient"
)

plt.title(
    "Position-specific association between transcript architecture and conservation"
)

plt.xticks(
    range(-20, 21, 2)
)

plt.legend()

plt.tight_layout()

plt.savefig(
    COMBINED_OUT,
    dpi=300
)

plt.close()

print("Saved:", COMBINED_OUT)
