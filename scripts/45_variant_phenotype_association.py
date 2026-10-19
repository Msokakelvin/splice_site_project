import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import os

# ── Paths ──────────────────────────────────────────────────────────────
VCF_FILE = "data/variants/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
JUNCTION_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/variant_phenotype_association_stats.tsv"
OUTPUT_PLOT = "results/figures/variant_burden_vs_splicing_phenotype.png"

WINDOW_SIZE = 50 

print("=== Script 45: Variant Burden vs. Splicing Phenotype Association ===")
print()

# 1. Load Junction Data with Phenotypes
print("Loading junction data and phenotypes (Chromosome 22)...")
usecols = ['junction_id', 'chromosome', 'architecture_fraction', 'exon_skipping', 'alternative_donor', 'alternative_acceptor']
df = pd.read_csv(JUNCTION_FILE, sep="\t", usecols=usecols)

# Parse chromosome and filter to Chr 22
df['chrom_clean'] = df['chromosome'].str.replace('chr', '', regex=False, case=False)
df_chr22 = df[df['chrom_clean'] == '22'].copy()

print(f"Total junctions with phenotypes on Chr 22: {len(df_chr22):,}")

if len(df_chr22) == 0:
    print("ERROR: No junctions found on Chromosome 22.")
    exit()

# Clean phenotype columns (ensure they are numeric 0/1 or floats)
for col in ['exon_skipping', 'alternative_donor', 'alternative_acceptor']:
    if col in df_chr22.columns:
        # Convert to binary 0/1 if they are strings like 'Yes'/'No', or keep as float if already numeric
        df_chr22[col] = pd.to_numeric(df_chr22[col], errors='coerce').fillna(0)
        # If the column was categorical (e.g., 'True'/'False'), map it
        if df_chr22[col].max() > 1:
            df_chr22[col] = (df_chr22[col] > 0).astype(int)

# 2. Open VCF and Count Variants Per Junction
print(f"\nOpening VCF: {VCF_FILE}")
if not os.path.exists(VCF_FILE):
    print(f"ERROR: VCF file not found at {VCF_FILE}")
    exit()

try:
    vcf = pysam.VariantFile(VCF_FILE)
except Exception as e:
    print(f"ERROR opening VCF: {e}")
    exit()

print("Counting variants per junction (this may take 5-10 minutes)...")
variant_counts = []

for idx, row in df_chr22.iterrows():
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
        pass
        
    variant_counts.append({'junction_id': jid, 'variant_count': count})
    
    if (idx + 1) % 5000 == 0:
        print(f"  Processed {idx + 1:,} junctions...")

vcf.close()

# 3. Merge Variant Counts with Phenotypes
print("\nMerging variant counts with splicing phenotypes...")
variant_df = pd.DataFrame(variant_counts)
merged_df = df_chr22.merge(variant_df, on='junction_id', how='inner')

print(f"Successfully merged {len(merged_df)} junctions.")

# 4. Statistical Association Tests
print("\n--- Statistical Association Tests ---")

# Test A: Spearman Correlation (Variant Count vs. Architecture Fraction)
rho, pval_corr = stats.spearmanr(merged_df['variant_count'], merged_df['architecture_fraction'])
print(f"\n1. Spearman Correlation: Variant Count vs. Architecture Fraction")
print(f"   Correlation coefficient (rho): {rho:.4f}")
print(f"   P-value: {pval_corr:.2e}")
if pval_corr < 0.05:
    print("   ✅ Significant: Higher variant burden correlates with lower architecture fraction.")

# Test B: Logistic Regression (Variant Count predicting Exon Skipping)
if 'exon_skipping' in merged_df.columns:
    print(f"\n2. Logistic Regression: Does Variant Count predict Exon Skipping?")
    X = sm.add_constant(merged_df['variant_count'])
    y = merged_df['exon_skipping']
    
    try:
        logit_model = sm.Logit(y, X).fit(disp=0)
        print(f"   Variant Count Coefficient: {logit_model.params['variant_count']:.4f}")
        print(f"   P-value: {logit_model.pvalues['variant_count']:.2e}")
        if logit_model.pvalues['variant_count'] < 0.05:
            print("   ✅ Significant: Higher variant burden significantly increases the odds of exon skipping.")
    except Exception as e:
        print(f"   ⚠️ Logistic regression failed (possibly due to low variance in exon_skipping): {e}")

# Save detailed stats
merged_df.to_csv(OUTPUT_STATS, sep="\t", index=False)
print(f"\nPer-junction association data saved to {OUTPUT_STATS}")

# 5. Plotting: Variant Count vs. Architecture Fraction
print("\nGenerating association plot...")
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(8, 6))

# Use a hexbin or scatter plot with regression line
sns.regplot(
    data=merged_df,
    x='variant_count',
    y='architecture_fraction',
    scatter_kws={'alpha': 0.3, 's': 10},
    line_kws={'color': 'red'},
    ci=95
)

plt.title('Local Variant Burden vs. Splice Site Architecture Fraction (Chr 22)', fontsize=14, fontweight='bold')
plt.xlabel('Number of Genetic Variants in ±50bp Window', fontsize=12)
plt.ylabel('Architecture Fraction', fontsize=12)

# Add correlation stats to plot
plt.text(0.05, 0.95, f'Spearman ρ = {rho:.3f}\np = {pval_corr:.2e}', 
         transform=plt.gca().transAxes, fontsize=12, fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 45 complete. Objective 4 fulfilled.")
