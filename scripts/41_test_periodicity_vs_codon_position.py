import pandas as pd
import numpy as np

INPUT_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_FILE = "data/processed/periodicity_analysis_summary.tsv"

print("=== Script 41: Testing 3-bp Periodicity vs Codon Position ===")
print()

print("Loading final table...")
df = pd.read_csv(INPUT_FILE, sep="\t")
print("  Loaded", len(df), "rows")

# Define codon dominance: which codon position has the highest count?
def get_dominant_codon(row):
    if row['cds'] == 0:
        return 'non_cds'
    codon_counts = [row['codon1'], row['codon2'], row['codon3']]
    max_idx = np.argmax(codon_counts)
    return f'codon{max_idx + 1}'

print("Assigning dominant codon position...")
df['dominant_codon'] = df.apply(get_dominant_codon, axis=1)

# Group by site, position, and dominant codon
print("Aggregating by site, position, and codon...")
summary = df.groupby(['site', 'position', 'dominant_codon']).agg({
    'architecture_fraction': ['mean', 'median', 'count'],
    'cds_frac': 'mean'
}).reset_index()

# Flatten multi-level columns
summary.columns = ['site', 'position', 'dominant_codon', 
                   'mean_arch_frac', 'median_arch_frac', 'n_junctions',
                   'mean_cds_frac']

print("Writing summary to", OUTPUT_FILE, "...")
summary.to_csv(OUTPUT_FILE, sep="\t", index=False)

print()
print("Done! Script 41 complete.")
print()
print("Summary statistics:")
print(summary.groupby('dominant_codon')['mean_arch_frac'].mean())
