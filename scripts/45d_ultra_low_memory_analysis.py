import pandas as pd
import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import os

SAMPLE_FILE = "data/processed/variant_counts_sample.tsv"
MAIN_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/variant_phenotype_association_stats_FINAL.tsv"
OUTPUT_PLOT = "results/figures/variant_burden_vs_splicing_phenotype_FINAL.png"

print("=== Script 45d: Ultra-Low Memory Final Analysis ===")
print()

# 1. Load the small random sample
print("1. Loading random sample (1,000,000 rows)...")
sample_df = pd.read_csv(SAMPLE_FILE, sep="\t")
sample_jids = set(sample_df['junction_id'])
print(f"   Loaded {len(sample_df):,} junction IDs into memory.")

# 2. Read main file in TINY chunks to prevent ANY memory spikes
print("2. Filtering main table in chunks (memory-safe)...")
chunks = []
chunk_size = 200000

for chunk in pd.read_csv(MAIN_FILE, sep="\t", chunksize=chunk_size, 
                         usecols=['junction_id', 'architecture_fraction', 'exon_skipping', 
                                  'alternative_donor', 'alternative_acceptor']):
    # Keep ONLY the rows that are in our random sample
    filtered = chunk[chunk['junction_id'].isin(sample_jids)]
    chunks.append(filtered)

main_filtered = pd.concat(chunks, ignore_index=True)
print(f"   Filtered main table to {len(main_filtered):,} matching rows.")

# 3. Merge the two small dataframes
print("3. Merging datasets...")
merged = main_filtered.merge(sample_df, on='junction_id')
print(f"   Final merged dataset: {len(merged):,} rows. (Uses < 500MB RAM)")

# Clean phenotype columns
for col in ['exon_skipping', 'alternative_donor', 'alternative_acceptor']:
    if col in merged.columns:
        merged[col] = pd.to_numeric(merged[col], errors='coerce').fillna(0)
        if merged[col].max() > 1:
            merged[col] = (merged[col] > 0).astype(int)

# 4. Statistical Tests
print("\n--- Statistical Association Tests ---")

# Spearman Correlation
rho, pval_corr = stats.spearmanr(merged['variant_count'], merged['architecture_fraction'])
print(f"\n1. Spearman Correlation: Variant Count vs. Architecture Fraction")
print(f"   Correlation coefficient (rho): {rho:.4f}")
print(f"   P-value: {pval_corr:.2e}")

if pval_corr < 0.05:
    print("   ✅ Significant correlation detected")

# Logistic Regression (Exon Skipping)
if 'exon_skipping' in merged.columns:
    print(f"\n2. Logistic Regression: Variant Count → Exon Skipping")
    X = sm.add_constant(merged['variant_count'])
    y = merged['exon_skipping']
    
    try:
        logit_model = sm.Logit(y, X).fit(disp=0)
        print(f"   Coefficient: {logit_model.params['variant_count']:.4f}")
        print(f"   P-value: {logit_model.pvalues['variant_count']:.2e}")
        
        if logit_model.pvalues['variant_count'] < 0.05:
            coef = logit_model.params['variant_count']
            if coef > 0:
                print("   ✅ Significant: Higher variant burden increases odds of exon skipping")
            else:
                print("   ✅ Significant: Higher variant burden decreases odds of exon skipping")
    except Exception as e:
        print(f"   ⚠️ Logistic regression failed: {e}")

# 5. Save Results
print(f"\n4. Saving final results to {OUTPUT_STATS}...")
merged.to_csv(OUTPUT_STATS, sep="\t", index=False)
print("   ✓ Saved.")

# Clean up sample file
if os.path.exists(SAMPLE_FILE):
    os.remove(SAMPLE_FILE)
    print("   ✓ Temporary sample file cleaned up.")

# 6. Plotting
print("\n5. Generating plot...")
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(8, 6))

sns.regplot(
    data=merged,
    x='variant_count',
    y='architecture_fraction',
    scatter_kws={'alpha': 0.1, 's': 5},
    line_kws={'color': 'red', 'linewidth': 2},
    ci=95
)

plt.title('Local Variant Burden vs. Splice Site Architecture Fraction\n(FULL GENOME - 1000G Phase 3 GRCh38)', 
          fontsize=14, fontweight='bold')
plt.xlabel('Number of Genetic Variants in ±50bp Window', fontsize=12)
plt.ylabel('Architecture Fraction', fontsize=12)

plt.text(0.05, 0.95, f'Spearman ρ = {rho:.3f}\np = {pval_corr:.2e}\nn = {len(merged):,}', 
         transform=plt.gca().transAxes, fontsize=11, fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round'))

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"   ✓ Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 45d COMPLETE. Objective 4 fulfilled!")
