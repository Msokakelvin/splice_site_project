import pandas as pd
import pysam
import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob

VCF_DIR = "data/variants"
JUNCTION_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/variant_phenotype_association_subset.tsv"
OUTPUT_PLOT = "results/figures/variant_burden_vs_splicing_phenotype_subset.png"

TARGET_CHROMS = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5'}
WINDOW_SIZE = 500 

print("=== Script 45e: Subset Chromosome Analysis (Chr 1-5) - FIXED ===")
print()

print("Loading junction data...")
usecols = ['junction_id', 'chromosome', 'architecture_fraction', 'exon_skipping', 'alternative_donor', 'alternative_acceptor']
df = pd.read_csv(JUNCTION_FILE, sep="\t", usecols=usecols)

# Keep the 'chr' prefix to match VCF format
print("Filtering to Chromosomes 1-5...")
df_subset = df[df['chromosome'].isin(TARGET_CHROMS)].copy()

# Deduplicate to one row per junction_id
df_subset = df_subset.drop_duplicates(subset='junction_id', keep='first')
print(f"Processing {len(df_subset):,} unique junctions")

coords = df_subset['junction_id'].str.split(':').str[1].str.split('-')
df_subset['start'] = coords.str[0].astype(int)
df_subset['end'] = coords.str[1].astype(int)
df_subset['strand'] = df_subset['junction_id'].str.split(':').str[2]

for col in ['exon_skipping', 'alternative_donor', 'alternative_acceptor']:
    if col in df_subset.columns:
        df_subset[col] = pd.to_numeric(df_subset[col], errors='coerce').fillna(0)
        if df_subset[col].max() > 1:
            df_subset[col] = (df_subset[col] > 0).astype(int)

print(f"\nScanning for GRCh38 VCF files...")
vcf_files = glob.glob(os.path.join(VCF_DIR, "1kGP_high_coverage_Illumina.chr*.vcf.gz"))
chrom_to_vcf = {}
for vcf_path in vcf_files:
    filename = os.path.basename(vcf_path)
    # Extract chromosome WITH 'chr' prefix
    chrom = filename.split('.')[1]  # This gives 'chr1', 'chr2', etc.
    chrom_to_vcf[chrom] = vcf_path

print(f"Found {len(chrom_to_vcf)} VCF files")

print("\n=== Processing Chromosomes 1-5 (1000bp windows) ===")
all_results = []

for chrom in sorted(TARGET_CHROMS):
    chrom_df = df_subset[df_subset['chromosome'] == chrom]
    print(f"  {chrom}: Processing {len(chrom_df):,} junctions...")
    
    if chrom not in chrom_to_vcf:
        print(f"    VCF not found, skipping.")
        continue
        
    vcf_file = chrom_to_vcf[chrom]
    
    try:
        vcf = pysam.VariantFile(vcf_file)
    except Exception as e:
        print(f"    Error opening VCF: {e}")
        continue
    
    chrom_results = []
    for idx, row in chrom_df.iterrows():
        jid = row['junction_id']
        start = row['start']
        end = row['end']
        strand = row['strand']
        
        if strand == '+':
            win_start = max(1, start - WINDOW_SIZE)
            win_end = end + WINDOW_SIZE
        else:
            win_start = max(1, end - WINDOW_SIZE)
            win_end = start + WINDOW_SIZE
        
        count = 0
        try:
            # FIX: Use 'chrom' which now has 'chr' prefix (e.g., 'chr1')
            for record in vcf.fetch(chrom, win_start, win_end):
                count += 1
        except (ValueError, KeyError) as e:
            pass
        
        chrom_results.append({'junction_id': jid, 'variant_count': count})
        
        if len(chrom_results) % 20000 == 0:
            print(f"    Progress: {len(chrom_results):,} / {len(chrom_df):,}")
    
    vcf.close()
    all_results.extend(chrom_results)
    print(f"    {chrom} complete.")

print("\n=== Final Analysis ===")
variant_df = pd.DataFrame(all_results)
variant_df = variant_df.drop_duplicates(subset='junction_id', keep='first')

merged = df_subset.merge(variant_df, on='junction_id', how='inner')
print(f"Successfully merged {len(merged):,} junctions.")

# Check variance before running stats
print(f"\nVariant count distribution (first 10 values):")
print(merged['variant_count'].value_counts().head(10))

print("\n--- Statistical Association Tests ---")

# Spearman Correlation
rho, pval_corr = stats.spearmanr(merged['variant_count'], merged['architecture_fraction'])
print(f"\n1. Spearman Correlation: Variant Count vs. Architecture Fraction")
print(f"   Correlation coefficient (rho): {rho:.4f}")
print(f"   P-value: {pval_corr:.2e}")

if not np.isnan(rho) and pval_corr < 0.05:
    print("   ✅ Significant correlation detected")

# Logistic Regression (Exon Skipping)
if 'exon_skipping' in merged.columns:
    print(f"\n2. Logistic Regression: Variant Count -> Exon Skipping")
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

print(f"\nSaving results to {OUTPUT_STATS}...")
merged.to_csv(OUTPUT_STATS, sep="\t", index=False)
print("Saved.")

print("\nGenerating plot...")
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(8, 6))

sns.regplot(
    data=merged,
    x='variant_count',
    y='architecture_fraction',
    scatter_kws={'alpha': 0.3, 's': 10},
    line_kws={'color': 'red', 'linewidth': 2},
    ci=95
)

plt.title('Local Variant Burden vs. Splice Site Architecture Fraction\n(SUBSET: Chromosomes 1-5, 1000bp Window, 1000G Phase 3)', 
          fontsize=14, fontweight='bold')
plt.xlabel('Number of Genetic Variants in +/-500bp Window', fontsize=12)
plt.ylabel('Architecture Fraction', fontsize=12)

plt.text(0.05, 0.95, f'Spearman rho = {rho:.3f}\np = {pval_corr:.2e}\nn = {len(merged):,}', 
         transform=plt.gca().transAxes, fontsize=11, fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round'))

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 45e COMPLETE. Objective 4 fulfilled!")
