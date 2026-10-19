import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import os
import glob

# ── Paths ──────────────────────────────────────────────────────────────
VCF_DIR = "data/variants"
JUNCTION_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/variant_phenotype_association_stats_FULL_GRCh38.tsv"
OUTPUT_PLOT = "results/figures/variant_burden_vs_splicing_phenotype_FULL_GRCh38.png"

WINDOW_SIZE = 50 

print("=== Script 45 FULL: Variant Burden vs. Splicing Phenotype (GRCh38 Full Genome) ===")
print()

# 1. Load Junction Data
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

# 2. Find all GRCh38 VCF files
print(f"\nScanning for GRCh38 VCF files in {VCF_DIR}...")
vcf_files = glob.glob(os.path.join(VCF_DIR, "1kGP_high_coverage_Illumina.chr*.vcf.gz"))
print(f"Found {len(vcf_files)} VCF files")

if len(vcf_files) == 0:
    print("ERROR: No GRCh38 VCF files found. Please run the download script first.")
    exit()

# Create mapping: chromosome -> VCF file
chrom_to_vcf = {}
for vcf_path in vcf_files:
    filename = os.path.basename(vcf_path)
    # Extract chromosome: 1kGP_high_coverage_Illumina.chr1.filtered... -> "1"
    chrom = filename.split('.')[2].replace('chr', '')
    chrom_to_vcf[chrom] = vcf_path

print(f"Successfully mapped {len(chrom_to_vcf)} chromosomes.")

# 3. Count Variants Per Junction
print("\nCounting variants per junction across the FULL GENOME...")
variant_counts = []

for idx, row in df.iterrows():
    jid = row['junction_id']
    chrom = row['chrom_clean']
    
    if chrom not in chrom_to_vcf:
        variant_counts.append({'junction_id': jid, 'variant_count': 0})
        continue
    
    vcf_file = chrom_to_vcf[chrom]
    
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
        with pysam.VariantFile(vcf_file) as vcf:
            for record in vcf.fetch(chrom, win_start, win_end):
                count += 1
    except (ValueError, KeyError, OSError):
        pass
        
    variant_counts.append({'junction_id': jid, 'variant_count': count})
    
    if (idx + 1) % 500000 == 0:
        print(f"  Processed {idx + 1:,} junctions ({(idx+1)/len(df)*100:.1f}%)...")

# 4. Merge and Analyze
print("\nMerging variant counts with splicing phenotypes...")
variant_df = pd.DataFrame(variant_counts)
merged_df = df.merge(variant_df, on='junction_id', how='inner')

print(f"Successfully merged {len(merged_df):,} junctions.")

print("\n--- Statistical Association Tests (Full Genome GRCh38) ---")

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
print(f"\nData saved to {OUTPUT_STATS}")

print("\nGenerating association plot...")
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(8, 6))

sns.regplot(
    data=merged_df,
    x='variant_count',
    y='architecture_fraction',
    scatter_kws={'alpha': 0.05, 's': 3},
    line_kws={'color': 'red', 'linewidth': 2},
    ci=95
)

plt.title('Local Variant Burden vs. Splice Site Architecture Fraction\n(FULL GENOME - 1000G Phase 3 GRCh38)', fontsize=14, fontweight='bold')
plt.xlabel('Number of Genetic Variants in ±50bp Window', fontsize=12)
plt.ylabel('Architecture Fraction', fontsize=12)

plt.text(0.05, 0.95, f'Spearman ρ = {rho:.3f}\np = {pval_corr:.2e}\nn = {len(merged_df):,}', 
         transform=plt.gca().transAxes, fontsize=11, fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round'))

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 45 FULL GRCh38 complete.")
