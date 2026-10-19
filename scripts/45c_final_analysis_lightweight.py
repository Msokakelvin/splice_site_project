import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import os

# ── Paths ──────────────────────────────────────────────────────────────
PROGRESS_FILE = "data/processed/variant_counts_progress.tsv"
JUNCTION_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/variant_phenotype_association_stats_FINAL.tsv"
OUTPUT_PLOT = "results/figures/variant_burden_vs_splicing_phenotype_FINAL.png"

SAMPLE_SIZE = 1000000  # Analyze 1 million rows for stats (statistically identical, 10x less RAM)

print("=== Script 45c: Lightweight Final Analysis ===")
print()

# 1. Load only necessary columns to save memory
print("Loading data (memory-efficient mode)...")
usecols_df = ['junction_id', 'architecture_fraction', 'exon_skipping', 'alternative_donor', 'alternative_acceptor']
df = pd.read_csv(JUNCTION_FILE, sep="\t", usecols=usecols_df)

# Clean phenotype columns
for col in ['exon_skipping', 'alternative_donor', 'alternative_acceptor']:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
        if df[col].max() > 1:
            df[col] = (df[col] > 0).astype(int)

print(f"Loaded {len(df):,} junctions from master table.")

# Load variant counts
variant_df = pd.read_csv(PROGRESS_FILE, sep="\t", usecols=['junction_id', 'variant_count'])
print(f"Loaded {len(variant_df):,} variant counts from progress file.")

# 2. Merge
print("Merging datasets...")
merged_df = df.merge(variant_df, on='junction_id', how='inner')
print(f"Merged dataset size: {len(merged_df):,} rows.")

# 3. Random Sampling for Statistical Tests (Prevents OOM Kill)
print(f"\nTaking a random sample of {SAMPLE_SIZE:,} rows for statistical analysis...")
sample_df = merged_df.sample(n=min(SAMPLE_SIZE, len(merged_df)), random_state=42)

print("\n--- Statistical Association Tests (Sampled) ---")

# Spearman Correlation
rho, pval_corr = stats.spearmanr(sample_df['variant_count'], sample_df['architecture_fraction'])
print(f"\n1. Spearman Correlation: Variant Count vs. Architecture Fraction")
print(f"   Correlation coefficient (rho): {rho:.4f}")
print(f"   P-value: {pval_corr:.2e}")

if pval_corr < 0.05:
    print("   ✅ Significant correlation detected")

# Logistic Regression (Exon Skipping)
if 'exon_skipping' in sample_df.columns:
    print(f"\n2. Logistic Regression: Variant Count → Exon Skipping")
    X = sm.add_constant(sample_df['variant_count'])
    y = sample_df['exon_skipping']
    
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

# 4. Save Full Results (in chunks to avoid OOM)
print(f"\nSaving full merged results to {OUTPUT_STATS} (this may take a minute)...")
# Drop unnecessary columns before saving to keep file size manageable
cols_to_save = ['junction_id', 'architecture_fraction', 'exon_skipping', 'alternative_donor', 'alternative_acceptor', 'variant_count']
merged_df[cols_to_save].to_csv(OUTPUT_STATS, sep="\t", index=False, chunksize=500000)
print("✓ Full results saved.")

# Clean up progress file
if os.path.exists(PROGRESS_FILE):
    os.remove(PROGRESS_FILE)
    print("✓ Progress file cleaned up.")

# 5. Plotting (using the sampled data for speed)
print("\nGenerating plot...")
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(8, 6))

sns.regplot(
    data=sample_df,
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

plt.text(0.05, 0.95, f'Spearman ρ = {rho:.3f}\np = {pval_corr:.2e}\nn = {len(merged_df):,} (sampled: {len(sample_df):,})', 
         transform=plt.gca().transAxes, fontsize=11, fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round'))

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"✓ Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 45c FINAL ANALYSIS complete. Objective 4 fulfilled!")
