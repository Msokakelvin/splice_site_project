import pandas as pd
import pysam
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# ── Paths ──────────────────────────────────────────────────────────────
VCF_FILE = "data/variants/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
JUNCTION_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_STATS = "data/processed/variant_density_stats_chr22.tsv"
OUTPUT_PLOT = "results/figures/variant_density_constitutive_vs_alternative_chr22.png"

# Window to check around each splice site (e.g., -50 to +50 bp)
WINDOW_SIZE = 50 

print("=== Script 44: Variant Density Analysis (Chromosome 22) ===")
print()

# 1. Load Junction Data
print("Loading junction data and filtering to Chromosome 22...")
df = pd.read_csv(JUNCTION_FILE, sep="\t", usecols=['junction_id', 'chromosome', 'architecture_fraction'])

# Parse chromosome to match VCF format (e.g., "chr22" -> "22")
df['chrom_clean'] = df['chromosome'].str.replace('chr', '', regex=False, case=False)

# Filter to Chromosome 22 only
df_chr22 = df[df['chrom_clean'] == '22'].copy()

print(f"Total junctions in dataset: {len(df):,}")
print(f"Junctions on Chromosome 22: {len(df_chr22):,}")

if len(df_chr22) == 0:
    print("ERROR: No junctions found on Chromosome 22. Check chromosome naming.")
    exit()

# Define Constitutive vs Alternative
df_chr22['junction_type'] = df_chr22['architecture_fraction'].apply(
    lambda x: 'Constitutive' if x >= 0.90 else 'Alternative'
)

# 2. Open VCF file
print(f"\nOpening VCF: {VCF_FILE}")
if not os.path.exists(VCF_FILE):
    print(f"ERROR: VCF file not found at {VCF_FILE}")
    print("Please ensure the Chr 22 VCF is downloaded and named correctly.")
    exit()

try:
    vcf = pysam.VariantFile(VCF_FILE)
except Exception as e:
    print(f"ERROR opening VCF: {e}")
    exit()

# 3. Count variants in windows
print("Counting variants in +/- 50bp windows...")
variant_counts = []

for idx, row in df_chr22.iterrows():
    jid = row['junction_id']
    jtype = row['junction_type']
    chrom = row['chrom_clean'] # This will be '22'
    
    # Parse junction coordinates: chr22:12345-67890:+
    parts = jid.split(':')
    coords = parts[1].split('-')
    start = int(coords[0])
    end = int(coords[1])
    strand = parts[2]
    
    # Define window: 50bp upstream of donor, 50bp downstream of acceptor
    if strand == '+':
        win_start = start - WINDOW_SIZE
        win_end = end + WINDOW_SIZE
    else:
        win_start = end - WINDOW_SIZE
        win_end = start + WINDOW_SIZE
        
    if win_start < 1:
        win_start = 1
        
    # Query VCF for this region
    count = 0
    try:
        for record in vcf.fetch(chrom, win_start, win_end):
            count += 1
    except ValueError:
        # Region might be out of bounds, count = 0
        pass
    except KeyError:
        # Chromosome not in VCF, count = 0
        pass
        
    variant_counts.append({
        'junction_id': jid,
        'junction_type': jtype,
        'variant_count': count,
        'window_size_bp': WINDOW_SIZE * 2
    })
    
    # Progress indicator
    if (idx + 1) % 5000 == 0:
        print(f"  Processed {idx + 1:,} junctions...")

vcf.close()

# 4. Analyze results
print("\nAnalyzing variant density...")
results_df = pd.DataFrame(variant_counts)

# Calculate density (variants per kb)
results_df['variant_density_per_kb'] = (results_df['variant_count'] / (results_df['window_size_bp'] / 1000))

# Statistical Test (Mann-Whitney U)
const_variants = results_df[results_df['junction_type'] == 'Constitutive']['variant_density_per_kb']
alt_variants = results_df[results_df['junction_type'] == 'Alternative']['variant_density_per_kb']

stat, p_value = stats.mannwhitneyu(const_variants, alt_variants, alternative='less')

print("\n--- Variant Density Summary (Chromosome 22) ---")
print(f"Constitutive mean variant density: {const_variants.mean():.2f} variants/kb")
print(f"Alternative mean variant density:  {alt_variants.mean():.2f} variants/kb")
print(f"P-value (Mann-Whitney U): {p_value:.2e}")

if p_value < 0.05:
    print("✅ Result: Constitutive junctions have significantly FEWER variants than alternative junctions.")
else:
    print("⚠️ Result: No significant difference in variant density.")

# Save stats
summary_stats = results_df.groupby('junction_type')['variant_density_per_kb'].agg(['mean', 'median', 'count', 'sum']).reset_index()
summary_stats.to_csv(OUTPUT_STATS, sep="\t", index=False)
print(f"\nStats saved to {OUTPUT_STATS}")

# 5. Plotting
print("\nGenerating plot...")
sns.set_theme(style="whitegrid", font_scale=1.1)
plt.figure(figsize=(8, 6))

sns.boxplot(
    data=results_df,
    x='junction_type',
    y='variant_density_per_kb',
    palette='Set2',
    order=['Constitutive', 'Alternative']
)

plt.title('Genetic Variant Density: Constitutive vs. Alternative Splice Sites (Chr 22)', fontsize=14, fontweight='bold')
plt.xlabel('Junction Type', fontsize=12)
plt.ylabel('Variant Density (variants per kb)', fontsize=12)

# Add p-value annotation
max_y = results_df['variant_density_per_kb'].max()
plt.text(0.5, max_y * 0.9, f'p = {p_value:.2e}', ha='center', va='center', fontsize=12, fontweight='bold', color='red')

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"Plot saved to {OUTPUT_PLOT}")
print("\n✅ Script 44 complete.")
