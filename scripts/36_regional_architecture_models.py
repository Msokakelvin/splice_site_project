import pandas as pd
import numpy as np
import statsmodels.formula.api as smf


INPUT = "data/processed/correct_phyloP_profiles.tsv"
MASTER = "data/processed/master_gene_junction_table.tsv"

OUTPUT = (
    "results/tables/"
    "regional_architecture_conservation.tsv"
)


print("Loading data...")

profiles = pd.read_csv(INPUT, sep="\t")
master = pd.read_csv(MASTER, sep="\t")

print("Profile rows:", len(profiles))
print("Master rows:", len(master))


# =========================================================
# 1. Restrict to multi-architecture genes
# =========================================================

master = master[
    master["total_architectures"] > 1
].copy()

print(
    "Multi-architecture rows:",
    len(master)
)

print(
    "Multi-architecture genes:",
    master["gene_id"].nunique()
)


# =========================================================
# 2. Intron length
# =========================================================

master["intron_length"] = (
    master["intron_end"]
    - master["intron_start"]
    + 1
)

master["log_intron_length"] = np.log10(
    master["intron_length"]
)


# =========================================================
# 3. Splice class
# =========================================================

master["splice_class_simple"] = np.where(
    master["splice_class"] == "canonical_GT_AG",
    "canonical",
    "noncanonical"
)


# =========================================================
# 4. Within-gene architecture fraction
# =========================================================

gene_mean = (
    master.groupby("gene_id")[
        "architecture_fraction"
    ].transform("mean")
)

master["within_gene_architecture"] = (
    master["architecture_fraction"]
    - gene_mean
)


# =========================================================
# 5. Merge profiles with model variables
# =========================================================

variables = master[
    [
        "gene_id",
        "junction_id",
        "within_gene_architecture",
        "log_intron_length",
        "splice_class_simple"
    ]
].copy()


df = profiles.merge(
    variables,
    on=["gene_id", "junction_id"],
    how="inner"
)

print("Merged rows:", len(df))
print("Merged genes:", df["gene_id"].nunique())


# =========================================================
# 6. Convert profiles to numerical lists
# =========================================================

def parse_profile(x):

    values = str(x).split(",")

    return [
        float(v)
        for v in values
        if v != ""
    ]


df["donor_profile"] = (
    df["donor_profile"].apply(parse_profile)
)

df["acceptor_profile"] = (
    df["acceptor_profile"].apply(parse_profile)
)


# =========================================================
# 7. Position mapping
# =========================================================

positions = (
    list(range(-20, 0))
    + list(range(1, 21))
)


# =========================================================
# 8. Calculate regional means
# =========================================================

def regional_mean(profile, selected_positions):

    if len(profile) != 40:
        return np.nan

    values = []

    for position in selected_positions:

        index = positions.index(position)

        values.append(profile[index])

    return np.mean(values)


# ---------------------------------------------------------
# DONOR REGIONS
# ---------------------------------------------------------

df["donor_exonic_flank"] = df[
    "donor_profile"
].apply(
    lambda x: regional_mean(
        x,
        list(range(-20, -2))
    )
)

df["donor_boundary"] = df[
    "donor_profile"
].apply(
    lambda x: regional_mean(
        x,
        [-2, -1, 1, 2]
    )
)

df["donor_intronic_flank"] = df[
    "donor_profile"
].apply(
    lambda x: regional_mean(
        x,
        list(range(3, 21))
    )
)


# ---------------------------------------------------------
# ACCEPTOR REGIONS
# ---------------------------------------------------------

df["acceptor_intronic_flank"] = df[
    "acceptor_profile"
].apply(
    lambda x: regional_mean(
        x,
        list(range(-20, -2))
    )
)

df["acceptor_boundary"] = df[
    "acceptor_profile"
].apply(
    lambda x: regional_mean(
        x,
        [-2, -1, 1, 2]
    )
)

df["acceptor_exonic_flank"] = df[
    "acceptor_profile"
].apply(
    lambda x: regional_mean(
        x,
        list(range(3, 21))
    )
)


# =========================================================
# 9. Run cluster-robust models
# =========================================================

results = []


def run_model(site, region, outcome):

    required = [
        outcome,
        "within_gene_architecture",
        "log_intron_length",
        "splice_class_simple",
        "gene_id"
    ]

    model_df = df[required].dropna().copy()

    print()
    print("--------------------------------------------")
    print("Site:", site)
    print("Region:", region)
    print("Outcome:", outcome)
    print("Rows:", len(model_df))
    print("Genes:", model_df["gene_id"].nunique())

    model = smf.ols(
        f"{outcome} ~ "
        "within_gene_architecture + "
        "log_intron_length + "
        "C(splice_class_simple)",
        data=model_df
    )

    result = model.fit(
        cov_type="cluster",
        cov_kwds={
            "groups": model_df["gene_id"]
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

    ci_low, ci_high = (
        result.conf_int()
        .loc["within_gene_architecture"]
    )

    results.append({
        "site": site,
        "region": region,
        "n": len(model_df),
        "genes": model_df["gene_id"].nunique(),
        "coefficient": coefficient,
        "SE": se,
        "CI_low": ci_low,
        "CI_high": ci_high,
        "p_value": p_value,
        "R_squared": result.rsquared
    })

    print(
        "Coefficient:",
        coefficient
    )

    print(
        "SE:",
        se
    )

    print(
        "95% CI:",
        ci_low,
        "-",
        ci_high
    )

    print(
        "P-value:",
        p_value
    )

    print(
        "R-squared:",
        result.rsquared
    )


# =========================================================
# 10. Donor models
# =========================================================

run_model(
    "donor",
    "exonic_flank",
    "donor_exonic_flank"
)

run_model(
    "donor",
    "boundary",
    "donor_boundary"
)

run_model(
    "donor",
    "intronic_flank",
    "donor_intronic_flank"
)


# =========================================================
# 11. Acceptor models
# =========================================================

run_model(
    "acceptor",
    "intronic_flank",
    "acceptor_intronic_flank"
)

run_model(
    "acceptor",
    "boundary",
    "acceptor_boundary"
)

run_model(
    "acceptor",
    "exonic_flank",
    "acceptor_exonic_flank"
)


# =========================================================
# 12. Save results
# =========================================================

results_df = pd.DataFrame(results)

results_df.to_csv(
    OUTPUT,
    sep="\t",
    index=False
)


# =========================================================
# 13. Final summary
# =========================================================

print()
print("============================================")
print("REGIONAL ANALYSIS COMPLETE")
print("============================================")

print()

print(
    results_df[
        [
            "site",
            "region",
            "coefficient",
            "CI_low",
            "CI_high",
            "p_value",
            "R_squared"
        ]
    ].to_string(index=False)
)

print()
print("Saved to:")
print(OUTPUT)
