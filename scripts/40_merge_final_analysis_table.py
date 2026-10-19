import pandas as pd

MASTER_FILE = "data/processed/master_gene_junction_table.tsv"
CDS_FILE    = "data/processed/junction_position_cds_summary.tsv"
OUTPUT_FILE = "data/processed/final_junction_analysis_table.tsv"

print("=== Script 40: Final Merge ===")

print("Loading master table...")
master_df = pd.read_csv(MASTER_FILE, sep="\t", dtype=str)
print("  Loaded", len(master_df), "junctions")

print("Loading CDS summary...")
cds_df = pd.read_csv(CDS_FILE, sep="\t")
print("  Loaded", len(cds_df), "rows")

print("Merging...")
final_df = pd.merge(master_df, cds_df, on="junction_id", how="inner")
print("  Retained", final_df["junction_id"].nunique(), "unique junctions")
print("  Total rows:", len(final_df))

print("Writing output...")
final_df.to_csv(OUTPUT_FILE, sep="\t", index=False)
print("Done!")
print("Columns:", list(final_df.columns))
