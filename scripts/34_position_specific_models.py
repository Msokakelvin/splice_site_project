import pandas as pd
import numpy as np
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

INPUT = "data/processed/correct_phyloP_profiles.tsv"
MASTER = "data/processed/master_gene_junction_table.tsv"

OUTPUT = "results/tables/position_specific_architecture_conservation.tsv"


print("Loading data...")

profiles = pd.read_csv(INPUT, sep="\t")
master = pd.read_csv(MASTER, sep="\t")

print("Profile rows:", len(profiles))
print("Master rows:", len(master))


# ---------------------------------------------------------
# 1. Keep only multi-architecture genes
# ---------------------------------------------------------

master = master[master["total_architectures"] > 1].copy()

print("Multi-architecture rows:", len(master))
print("Multi-architecture genes:", master["gene_id"].nunique())


# ---------------------------------------------------------
# 2. Calculate intron length
# ---------------------------------------------------------

master["intron_length"] = (
    master["intron_end"] - master["intron_start"] + 1
)

master["log_intron_length"] = np.log10(
    master["intron_length"]
)


# ---------------------------------------------------------
# 3. Simplify splice class
# ---------------------------------------------------------

master["splice_class_simple"] = np.where(
    master["splice_class"] == "canonical_GT_AG",
    "canonical",
    "noncanonical"
)


# ---------------------------------------------------------
# 4. Calculate within-gene architecture fraction
# ---------------------------------------------------------

gene_mean = (
    master.groupby("gene_id")["architecture_fraction"]
    .transform("mean")
)

master["within_gene_architecture"] = (
    master["architecture_fraction"] - gene_mean
)


# ---------------------------------------------------------
# 5. Merge profiles with model variables
# ---------------------------------------------------------

columns_needed = [
    "gene_id",
    "junction_id",
    "within_gene_architecture",
    "log_intron_length",
    "splice_class_simple"
]

model_variables = master[columns_needed].copy()

df = profiles.merge(
    model_variables,
    on=["gene_id", "junction_id"],
    how="inner"
)

print("Merged rows:", len(df))
print("Merged genes:", df["gene_id"].nunique())


# ---------------------------------------------------------
# 6. Convert profile strings into numeric arrays
# ---------------------------------------------------------

def parse_profile(x):

    values = str(x).split(",")

    return [
        float(v)
        for v in values
        if v != ""
    ]


df["donor_profile"] = df["donor_profile"].apply(parse_profile)
df["acceptor_profile"] = df["acceptor_profile"].apply(parse_profile)


# ---------------------------------------------------------
# 7. Define positions
# ---------------------------------------------------------

positions = list(range(-20, 0)) + list(range(1, 21))


# ---------------------------------------------------------
# 8. Run position-specific models
# ---------------------------------------------------------

results = []

def run_position(site, position, values):

    temp = df[
        [
            "gene_id",
            "within_gene_architecture",
            "log_intron_length",
            "splice_class_simple"
        ]
    ].copy()

    temp["phyloP"] = values

    required = [
        "phyloP",
        "within_gene_architecture",
        "log_intron_length",
        "splice_class_simple",
        "gene_id"
    ]

    temp = temp[required].dropna()

    if len(temp) < 1000:
        print(
            f"Skipping {site} {position}: "
            f"only {len(temp)} observations"
        )
        return

    model = smf.ols(
        "phyloP ~ within_gene_architecture "
        "+ log_intron_length "
        "+ C(splice_class_simple)",
        data=temp
    )

    result = model.fit(
        cov_type="cluster",
        cov_kwds={
            "groups": temp["gene_id"]
        }
    )

    coefficient = result.params[
        "within_gene_architecture"
    ]

    se = result.bse[
        "within_gene_architecture"
    ]

    p_value = result.pvalues[
        "within_gene_architecture"
    ]

    ci_low, ci_high = result.conf_int().loc[
        "within_gene_architecture"
    ]

    results.append({
        "site": site,
        "position": position,
        "n": len(temp),
        "genes": temp["gene_id"].nunique(),
        "coefficient": coefficient,
        "SE": se,
        "CI_low": ci_low,
        "CI_high": ci_high,
        "p_value": p_value,
        "R_squared": result.rsquared
    })


# ---------------------------------------------------------
# 9. Donor positions
# ---------------------------------------------------------

print("\nAnalyzing donor positions...")

for i, position in enumerate(positions):

    values = df["donor_profile"].apply(
        lambda x: x[i] if len(x) == 40 else np.nan
    )

    run_position(
        "donor",
        position,
        values
    )


# ---------------------------------------------------------
# 10. Acceptor positions
# ---------------------------------------------------------

print("\nAnalyzing acceptor positions...")

for i, position in enumerate(positions):

    values = df["acceptor_profile"].apply(
        lambda x: x[i] if len(x) == 40 else np.nan
    )

    run_position(
        "acceptor",
        position,
        values
    )


# ---------------------------------------------------------
# 11. Convert results to dataframe
# ---------------------------------------------------------

results_df = pd.DataFrame(results)


# ---------------------------------------------------------
# 12. Multiple-testing correction
# ---------------------------------------------------------

reject, fdr, _, _ = multipletests(
    results_df["p_value"],
    alpha=0.05,
    method="fdr_bh"
)

results_df["FDR"] = fdr
results_df["FDR_significant"] = reject


# ---------------------------------------------------------
# 13. Sort by site and position
# ---------------------------------------------------------

site_order = {
    "donor": 0,
    "acceptor": 1
}

results_df["site_order"] = (
    results_df["site"].map(site_order)
)

results_df = results_df.sort_values(
    ["site_order", "position"]
)

results_df = results_df.drop(
    columns=["site_order"]
)


# ---------------------------------------------------------
# 14. Save results
# ---------------------------------------------------------

results_df.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


# ---------------------------------------------------------
# 15. Print summary
# ---------------------------------------------------------

print("\n==============================================")
print("POSITION-SPECIFIC ANALYSIS COMPLETE")
print("==============================================")

print("Results:", len(results_df))

print(
    "\nFDR-significant positions:",
    results_df["FDR_significant"].sum()
)

print("\nSignificant positions:")

print(
    results_df[
        results_df["FDR_significant"]
    ][
        [
            "site",
            "position",
            "coefficient",
            "CI_low",
            "CI_high",
            "FDR"
        ]
    ].to_string(index=False)
)

print("\nSaved to:")
print(OUTPUT)
