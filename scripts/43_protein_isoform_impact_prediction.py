import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict

INPUT_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/protein_isoform_impact_stats.tsv"
OUTPUT_PLOT = "results/figures/protein_isoform_impact_distribution.png"

print("=== Script 43: Predicting Protein Isoform Impact (Memory-Efficient) ===")
print()

# 1. Load only the columns we need
print("Loading final analysis table (selective columns)...")
usecols = ['junction_id', 'architecture_fraction', 'cds_frac', 'cds', 'codon1', 'codon2', 'codon3']
df = pd.read_csv(INPUT_FILE, sep="\t", usecols=usecols)

print(f"Loaded {len(df):,} rows")
print()

# 2. Define Alternative Junctions (architecture_fraction < 0.90)
alt_mask = df['architecture_fraction'] < 0.90
alt_df = df[alt_mask].copy()

print(f"Total junction-position rows: {len(df):,}")
print(f"  - Constitutive (arch_frac >= 0.90): {(~alt_mask).sum():,}")
print(f"  - Alternative (arch_frac < 0.90): {len(alt_df):,}")
print()

# 3. Vectorized prediction (much faster than apply)
print("Predicting isoform impact (vectorized)...")

# Initialize impact column
alt_df['predicted_impact'] = "UTR / Non-coding Impact"

# Condition 1: If cds_frac >= 0.20, check codon phase
cds_mask = alt_df['cds_frac'] >= 0.20
cds_mask &= alt_df['cds'] > 0

# For CDS rows, calculate dominance ratio
codon_max = alt_df[['codon1', 'codon2', 'codon3']].max(axis=1)
dominance_ratio = codon_max / alt_df['cds'].replace(0, 1)  # avoid division by zero

# Condition 2: If dominance >= 0.60, it's in-frame
inframe_mask = cds_mask & (dominance_ratio >= 0.60)
alt_df.loc[inframe_mask, 'predicted_impact'] = "Potential In-Frame Isoform"

# Condition 3: If CDS but mixed codon phases, it's frameshift/NMD risk
frameshift_mask = cds_mask & (dominance_ratio < 0.60)
alt_df.loc[frameshift_mask, 'predicted_impact'] = "Frameshift / NMD Risk"

# 4. Aggregate statistics (group by junction_id to get unique counts)
print("Aggregating statistics...")

# Group by junction_id and predicted_impact
junction_impact = alt_df.groupby(['junction_id', 'predicted_impact']).size().reset_index(name='count')

# For each junction, take the most common impact (mode)
junction_mode = junction_impact.sort_values('count', ascending=False).drop_duplicates('junction_id')

# Count unique junctions per impact category
impact_counts = junction_mode['predicted_impact'].value_counts().reset_index()
impact_counts.columns = ['predicted_impact', 'unique_junctions']

# Calculate percentages
total_junctions = impact_counts['unique_junctions'].sum()
impact_counts['percentage'] = (impact_counts['unique_junctions'] / total_junctions * 100).round(1)

# Also calculate mean cds_frac for each category
mean_cds = alt_df.groupby('predicted_impact')['cds_frac'].mean().reset_index()
mean_cds.columns = ['predicted_impact', 'mean_cds_frac']

impact_summary = impact_counts.merge(mean_cds, on='predicted_impact')

print("\n--- Protein Isoform Impact Summary (Alternative Junctions) ---")
print(impact_summary.to_string(index=False))

impact_summary.to_csv(OUTPUT_STATS, sep="\t", index=False)
print(f"\nStats saved to {OUTPUT_STATS}")

# 5. Plotting
print("\nGenerating plot...")
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(10, 6))

# Order the categories logically
order = ["UTR / Non-coding Impact", "Potential In-Frame Isoform", "Frameshift / NMD Risk"]
plot_df = impact_summary.set_index('predicted_impact').reindex(order).reset_index()

sns.barplot(
    data=plot_df,
    x='predicted_impact',
    y='unique_junctions',
    palette=['#4C72B0', '#55A868', '#C44E52']
)

plt.title('Predicted Protein-Level Impact of Alternative Splice Sites', fontsize=14, fontweight='bold')
plt.xlabel('Predicted Impact Category', fontsize=12)
plt.ylabel('Number of Unique Alternative Junctions', fontsize=12)

# Add percentage labels on top of bars
for i, row in plot_df.iterrows():
    plt.text(i, row['unique_junctions'] + (plot_df['unique_junctions'].max() * 0.02), 
             f"{row['percentage']}%", 
             ha='center', va='bottom', fontsize=11, fontweight='bold')

plt.xticks(rotation=15, fontsize=11)
plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 43 complete.")
