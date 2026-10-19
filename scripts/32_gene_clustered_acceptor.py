import pandas as pd
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
# 4. Fit acceptor conservation model
# ---------------------------------------------------------

print("\nFitting gene-clustered acceptor model...")

model = smf.ols(
    "acceptor_motif_mean ~ within_gene_architecture + "
    "log_intron_length + C(splice_class_simple)",
    data=df
)

result = model.fit(
    cov_type="cluster",
    cov_kwds={"groups": df["gene_id"]}
)


# ---------------------------------------------------------
# 5. Display results
# ---------------------------------------------------------

print("\n==============================================")
print("GENE-CLUSTERED ACCEPTOR CONSERVATION MODEL")
print("==============================================")

print(result.summary())
