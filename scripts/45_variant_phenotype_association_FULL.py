import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import os

# ── Paths ──────────────────────────────────────────────────────────────
VCF_FILE = "data/variants/ALL.wgs.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
JUNCTION_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/variant_phenotype_association_stats_FULL.tsv"
OUTPUT_PLOT = "results/figures/variant_burden_vs_splicing_phenotype_FULL.png"

WINDOW_SIZE = 50 

print("=== Script 45 FULL: Variant Burden vs. Splicing Phenotype (Full Genome) ===")
print()

# 1. Load Junction Data with Phenotypes
print("Loading full genome junction data and phenotypes...")
usecols = ['junction_id', 'chromosome', 'architecture_fraction', 'exon_skipping', 'alternative_donor', 'alternative_acceptor']
df = pd.read_csv(JUNCTION_FILE, sep="\t", usecols=usecols)

# Parse chromosome (e.g., "chr1" -> "1", "chrX" -> "X")
df['chrom_clean'] = df['chromosome'].str.replace('chr', '', regex=False, case=False)

print(f"Total junctions to process: {len(df):,}")

# Clean phenotype columns
for col in ['exon_skipping', 'alternative_donor', 'alternative_acceptor']:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
        if df[col].max() > 1:
            df[col] = (df[col] > 0).astype(int)

# 2. Open VCF
print(f"\nOpening VCF: {VCF_FILE}")
if not os.path.exists(VCF_FILE):
    print(f"ERROR: Full genome VCF not found at {VCF_FILE}")
    exit()

try:
    vcf = pysam.VariantFile(VCF_FILE)
except Exception as e:
    print(f"ERROR opening VCF: {e}")
    exit()

# 3. Count Variants Per Junction (Full Genome)
print("Counting variants per junction across the FULL GENOME (this will take 1-2 hours)...")
variant_counts = []

for idx, row in df.iterrows():
    jid = row['junction_id']
    chrom = row['chrom_clean']
    
    parts = jid.split(':')
    coords = parts[1].split('-')
    start = int(coords[0])
    end = int(coords[1])
    strand = parts[2]
    
    if strand == '+':
        win_start = start - WINDOW_SIZE
        win_end = end + WINDOW_SIZE
    else:
        win_start = end - WINDOW_SIZE
        win_end = start + WINDOW_SIZE
        
    if win_start < 1:
        win_start = 1
        
    count = 0
    try:
        for record in vcf.fetch(chrom, win_start, win_end):
            count += 1
    except (ValueError, KeyError):
        # Out of bounds or unplaced scaffold
        pass
        
    variant_counts.append({'junction_id': jid, 'variant_count': count})
    
    # Progress indicator (every 200,000 junctions)
    if (idx + 1) % 200000 == 0:
        print(f"  Processed {idx + 1:,} junctions ({(idx+1)/len(df)*100:.1f}%)...")

vcf.close()

# 4. Merge and Analyze
print("\nMerging variant counts with splicing phenotypes...")
variant_df = pd.DataFrame(variant_counts)
merged_df = df.merge(variant_df, on='junction_id', how='inner')

print(f"Successfully merged {len(merged_df)} junctions.")

print("\n--- Statistical Association Tests (Full Genome) ---")

# Spearman Correlation
rho, pval_corr = stats.spearmanr(merged_df['variant_count'], merged_df['architecture_fraction'])
print(f"\n1. Spearman Correlation: Variant Count vs. Architecture Fraction")
print(f"   Correlation coefficient (rho): {rho:.4f}")
print(f"   P-value: {pval_corr:.2e}")

# Logistic Regression (Exon Skipping)
if 'exon_skipping' in merged_df.columns:
    print(f"\n2. Logistic Regression: Does Variant Count predict Exon Skipping?")
    X = sm.add_constant(merged_df['variant_count'])
    y = merged_df['exon_skipping']
    
    try:
        logit_model = sm.Logit(y, X).fit(disp=0)
        print(f"   Variant Count Coefficient: {logit_model.params['variant_count']:.4f}")
        print(f"   P-value: {logit_model.pvalues['variant_count']:.2e}")
    except Exception as e:
        print(f"   ⚠️ Logistic regression failed: {e}")

# Save and Plot
merged_df.to_csv(OUTPUT_STATS, sep="\t", index=False)
print(f"\nPer-junction association data saved to {OUTPUT_STATS}")

print("\nGenerating association plot...")
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(8, 6))

sns.regplot(
    data=merged_df,
    x='variant_count',
    y='architecture_fraction',
    scatter_kws={'alpha': 0.1, 's': 5}, # Lower alpha for 9.9M points
    line_kws={'color': 'red'},
    ci=95
)

plt.title('Local Variant Burden vs. Splice Site Architecture Fraction (FULL GENOME)', fontsize=14, fontweight='bold')
plt.xlabel('Number of Genetic Variants in ±50bp Window', fontsize=12)
plt.ylabel('Architecture Fraction', fontsize=12)

plt.text(0.05, 0.95, f'Spearman ρ = {rho:.3f}\np = {pval_corr:.2e}', 
         transform=plt.gca().transAxes, fontsize=12, fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 45 FULL complete.")
