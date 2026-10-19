import pandas as pd
import numpy as np
import statsmodels.formula.api as smf


# ---------------------------------------------------------
# 1. Input
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
# 3. Calculate within-gene architecture fraction
# ---------------------------------------------------------

gene_mean = (
    df.groupby("gene_id")["architecture_fraction"]
      .transform("mean")
)

df["within_gene_architecture"] = (
    df["architecture_fraction"] - gene_mean
)


# ---------------------------------------------------------
# 4. Fit gene-clustered model
# ---------------------------------------------------------

print("\nFitting gene-clustered donor model...")

model = smf.ols(
    "donor_motif_mean ~ within_gene_architecture + "
    "log_intron_length + C(splice_class_simple)",
    data=df
)

result = model.fit(
    cov_type="cluster",
    cov_kwds={"groups": df["gene_id"]}
)


# ---------------------------------------------------------
# 5. Results
# ---------------------------------------------------------

print("\n==============================================")
print("GENE-CLUSTERED DONOR CONSERVATION MODEL")
print("==============================================")

print(result.summary())
