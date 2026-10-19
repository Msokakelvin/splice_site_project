import pandas as pd
import statsmodels.formula.api as smf


# ---------------------------------------------------------
# 1. Input file
# ---------------------------------------------------------

INPUT = "data/processed/model_data.tsv"


# ---------------------------------------------------------
# 2. Load data
# ---------------------------------------------------------

print("Loading model data...")

df = pd.read_csv(INPUT, sep="\t")

print("Rows:", len(df))
print("Genes:", df["gene_id"].nunique())


# ---------------------------------------------------------
# 3. Fit mixed-effects model
# ---------------------------------------------------------

print("\nFitting donor conservation model...")

model = smf.mixedlm(
    "donor_motif_mean ~ architecture_fraction + log_intron_length + C(splice_class_simple)",
    data=df,
    groups=df["gene_id"]
)

result = model.fit(
    method="lbfgs",
    reml=False
)


# ---------------------------------------------------------
# 4. Display results
# ---------------------------------------------------------

print("\n==============================================")
print("DONOR CONSERVATION MIXED-EFFECTS MODEL")
print("==============================================")

print(result.summary())
