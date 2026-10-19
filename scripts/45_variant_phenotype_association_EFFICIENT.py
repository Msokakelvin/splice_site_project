import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import os
import glob
import gc
from collections import defaultdict

# ── Paths ──────────────────────────────────────────────────────────────
VCF_DIR = "data/variants"
JUNCTION_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/variant_phenotype_association_stats_EFFICIENT.tsv"
OUTPUT_PLOT = "results/figures/variant_burden_vs_splicing_phenotype_EFFICIENT.png"
PROGRESS_FILE = "data/processed/variant_counts_progress.tsv"  # Resume support

WINDOW_SIZE = 50 

print("=== Script 45 EFFICIENT: Chromosome-Grouped Analysis ===")
print()

# 1. Load Junction Data
print("Loading junction data...")
usecols = ['junction_id', 'chromosome', 'architecture_fraction', 'exon_skipping', 
           'alternative_donor', 'alternative_acceptor']
df = pd.read_csv(JUNCTION_FILE, sep="\t", usecols=usecols)
df['chrom_clean'] = df['chromosome'].str.replace('chr', '', regex=False, case=False)

# Parse coordinates once
print("Parsing junction coordinates...")
df['jid_parts'] = df['junction_id'].str.split(':')
df['chrom_from_jid'] = df['jid_parts'].str[0].str.replace('chr', '', regex=False)
coords = df['jid_parts'].str[1].str.split('-')
df['start'] = coords.str[0].astype(int)
df['end'] = coords.str[1].astype(int)
df['strand'] = df['jid_parts'].str[2]

# Clean phenotype columns
for col in ['exon_skipping', 'alternative_donor', 'alternative_acceptor']:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)
        if df[col].max() > 1:
            df[col] = (df[col] > 0).astype(int)

print(f"Total junctions: {len(df):,}")

# 2. Find VCF files
print(f"\nScanning for GRCh38 VCF files...")
vcf_files = glob.glob(os.path.join(VCF_DIR, "1kGP_high_coverage_Illumina.chr*.vcf.gz"))
print(f"Found {len(vcf_files)} VCF files")

chrom_to_vcf = {}
for vcf_path in vcf_files:
    filename = os.path.basename(vcf_path)
    chrom = filename.split('.')[1].replace('chr', '')
    chrom_to_vcf[chrom] = vcf_path

print(f"Mapped {len(chrom_to_vcf)} chromosomes")

# 3. Check for existing progress (resume support)
completed_junctions = set()
if os.path.exists(PROGRESS_FILE):
    print(f"\nFound existing progress file. Resuming...")
    progress_df = pd.read_csv(PROGRESS_FILE, sep="\t", usecols=['junction_id', 'variant_count'])
    completed_junctions = set(progress_df['junction_id'].values)
    all_results = progress_df.to_dict('records')
    print(f"  Already processed: {len(completed_junctions):,} junctions")
else:
    all_results = []

# 4. Process by Chromosome (EFFICIENT - opens VCF once per chromosome)
print("\n=== Processing by Chromosome ===")
chromosomes_to_process = sorted(chrom_to_vcf.keys(), 
                                 key=lambda x: (x.isdigit(), int(x) if x.isdigit() else x))

total_junctions = len(df)
total_processed = len(completed_junctions)

for chrom in chromosomes_to_process:
    # Filter junctions for this chromosome
    chrom_df = df[(df['chrom_clean'] == chrom) | (df['chrom_from_jid'] == chrom)].copy()
    
    # Skip already completed
    chrom_df = chrom_df[~chrom_df['junction_id'].isin(completed_junctions)]
    
    if len(chrom_df) == 0:
        print(f"  Chr {chrom}: Already complete, skipping")
        continue
    
    print(f"\n  Chr {chrom}: Processing {len(chrom_df):,} junctions...")
    
    vcf_file = chrom_to_vcf[chrom]
    
    # Open VCF ONCE for this chromosome
    try:
        vcf = pysam.VariantFile(vcf_file)
    except Exception as e:
        print(f"    ✗ Error opening VCF: {e}")
        continue
    
    # Process all junctions for this chromosome
    chrom_results = []
    processed_in_chrom = 0
    
    for idx, row in chrom_df.iterrows():
        jid = row['junction_id']
        start = row['start']
        end = row['end']
        strand = row['strand']
        
        # Calculate window
        if strand == '+':
            win_start = max(1, start - WINDOW_SIZE)
            win_end = end + WINDOW_SIZE
        else:
            win_start = max(1, end - WINDOW_SIZE)
            win_end = start + WINDOW_SIZE
        
        # Count variants
        count = 0
        try:
            for record in vcf.fetch(chrom, win_start, win_end):
                count += 1
        except (ValueError, KeyError):
            pass
        
        chrom_results.append({'junction_id': jid, 'variant_count': count})
        processed_in_chrom += 1
        
        # Progress update every 50,000 junctions
        if processed_in_chrom % 50000 == 0:
            total_processed += 50000
            print(f"    Progress: {processed_in_chrom:,} / {len(chrom_df):,} "
                  f"(Total: {total_processed:,} / {total_junctions:,})")
    
    vcf.close()
    
    # Save progress for this chromosome
    all_results.extend(chrom_results)
    progress_df = pd.DataFrame(all_results)
    progress_df.to_csv(PROGRESS_FILE, sep="\t", index=False)
    
    total_processed += len(chrom_df)
    print(f"    ✓ Chr {chrom} complete. Total: {total_processed:,} / {total_junctions:,}")
    
    # Free memory
    del chrom_df, chrom_results
    gc.collect()

# 5. Final Merge and Analysis
print("\n=== Final Analysis ===")
variant_df = pd.read_csv(PROGRESS_FILE, sep="\t")
merged_df = df.merge(variant_df, on='junction_id', how='inner')

print(f"Successfully merged {len(merged_df):,} junctions.")

print("\n--- Statistical Association Tests ---")

# Spearman Correlation
rho, pval_corr = stats.spearmanr(merged_df['variant_count'], merged_df['architecture_fraction'])
print(f"\n1. Spearman Correlation: Variant Count vs. Architecture Fraction")
print(f"   Correlation coefficient (rho): {rho:.4f}")
print(f"   P-value: {pval_corr:.2e}")

# Logistic Regression (Exon Skipping)
if 'exon_skipping' in merged_df.columns:
    print(f"\n2. Logistic Regression: Variant Count → Exon Skipping")
    X = sm.add_constant(merged_df['variant_count'])
    y = merged_df['exon_skipping']
    
    try:
        logit_model = sm.Logit(y, X).fit(disp=0)
        print(f"   Coefficient: {logit_model.params['variant_count']:.4f}")
        print(f"   P-value: {logit_model.pvalues['variant_count']:.2e}")
    except Exception as e:
        print(f"   ⚠️ Failed: {e}")

# Save final results
merged_df.to_csv(OUTPUT_STATS, sep="\t", index=False)
print(f"\n✓ Stats saved to {OUTPUT_STATS}")

# Clean up progress file
if os.path.exists(PROGRESS_FILE):
    os.remove(PROGRESS_FILE)
    print("✓ Progress file cleaned up")

# 6. Plotting
print("\nGenerating plot...")
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

plt.title('Local Variant Burden vs. Splice Site Architecture Fraction\n(FULL GENOME - 1000G Phase 3 GRCh38)', 
          fontsize=14, fontweight='bold')
plt.xlabel('Number of Genetic Variants in ±50bp Window', fontsize=12)
plt.ylabel('Architecture Fraction', fontsize=12)

plt.text(0.05, 0.95, f'Spearman ρ = {rho:.3f}\np = {pval_corr:.2e}\nn = {len(merged_df):,}', 
         transform=plt.gca().transAxes, fontsize=11, fontweight='bold',
         bbox=dict(facecolor='white', alpha=0.9, edgecolor='gray', boxstyle='round'))

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"✓ Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 45 EFFICIENT complete.")
