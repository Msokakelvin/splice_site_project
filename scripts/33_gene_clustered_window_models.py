import pandas as pd
import numpy as np
import statsmodels.formula.api as smf


# ---------------------------------------------------------
# 1. Input
# ---------------------------------------------------------

INPUT = "data/processed/architecture_conservation.tsv"


# ---------------------------------------------------------
# 2. Load data
# ---------------------------------------------------------

print("Loading architecture conservation data...")

df = pd.read_csv(INPUT, sep="\t")

# Keep only multi-architecture genes
df = df[df["total_architectures"] > 1].copy()

print("Rows:", len(df))
print("Genes:", df["gene_id"].nunique())


# ---------------------------------------------------------
# 3. Calculate intron length
# ---------------------------------------------------------

df["intron_length"] = (
    df["intron_end"] -
    df["intron_start"] +
    1
)

df["log_intron_length"] = np.log10(
    df["intron_length"]
)


# ---------------------------------------------------------
# 4. Simplify splice class
# ---------------------------------------------------------

df["splice_class_simple"] = np.where(
    df["splice_class"] == "canonical_GT_AG",
    "canonical",
    "noncanonical"
)


# ---------------------------------------------------------
# 5. Calculate within-gene architecture fraction
# ---------------------------------------------------------

gene_mean = (
    df.groupby("gene_id")["architecture_fraction"]
      .transform("mean")
)

df["within_gene_architecture"] = (
    df["architecture_fraction"] -
    gene_mean
)


# ---------------------------------------------------------
# 6. Function to fit clustered model
# ---------------------------------------------------------

def run_model(outcome):

    print("\n----------------------------------------------")
    print("Outcome:", outcome)
    print("----------------------------------------------")

    # Variables required for this model
    required = [
        outcome,
        "within_gene_architecture",
        "log_intron_length",
        "splice_class_simple",
        "gene_id"
    ]

    # Remove missing observations BEFORE fitting
    model_df = df[required].dropna().copy()

    print("Model rows:", len(model_df))
    print("Model genes:", model_df["gene_id"].nunique())

    # Fit ordinary least squares
    model = smf.ols(
        f"{outcome} ~ within_gene_architecture + "
        "log_intron_length + C(splice_class_simple)",
        data=model_df
    )

    # Cluster standard errors by gene
    result = model.fit(
        cov_type="cluster",
        cov_kwds={
            "groups": model_df["gene_id"]
        }
    )

    # Extract architecture effect
    coefficient = result.params[
        "within_gene_architecture"
    ]

    se = result.bse[
        "within_gene_architecture"
    ]

    p = result.pvalues[
        "within_gene_architecture"
    ]

    ci_low, ci_high = (
        result.conf_int()
        .loc["within_gene_architecture"]
    )

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
        p
    )

    print(
        "R-squared:",
        result.rsquared
    )

    return result


# ---------------------------------------------------------
# 7. Donor ±20 bp
# ---------------------------------------------------------

donor_result = run_model(
    "donor_window_mean"
)


# ---------------------------------------------------------
# 8. Acceptor ±20 bp
# ---------------------------------------------------------

acceptor_result = run_model(
    "acceptor_window_mean"
)


# ---------------------------------------------------------
# 9. Donor boundary
# ---------------------------------------------------------

donor_boundary_result = run_model(
    "donor_boundary_mean"
)


# ---------------------------------------------------------
# 10. Acceptor boundary
# ---------------------------------------------------------

acceptor_boundary_result = run_model(
    "acceptor_boundary_mean"
)
