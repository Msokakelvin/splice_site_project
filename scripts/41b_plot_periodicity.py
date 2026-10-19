import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

INPUT_FILE = "data/processed/periodicity_analysis_summary.tsv"
OUTPUT_PLOT = "results/periodicity_by_codon_position.png"

print("=== Script 41b: Plotting 3-bp Periodicity ===")

# Load data
df = pd.read_csv(INPUT_FILE, sep="\t")

# Ensure 'dominant_codon' is a categorical variable with a logical order
codon_order = ['non_cds', 'codon1', 'codon2', 'codon3']
df['dominant_codon'] = pd.Categorical(df['dominant_codon'], categories=codon_order, ordered=True)

# Set up the plotting style
sns.set_theme(style="whitegrid", font_scale=1.1)
fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

# Plot 1: Donor Site (-20 to -1)
donor_df = df[df['site'] == 'donor'].copy()
donor_df = donor_df.sort_values('position')

sns.lineplot(
    data=donor_df,
    x='position',
    y='mean_arch_frac',
    hue='dominant_codon',
    hue_order=codon_order,
    marker='o',
    ax=axes[0],
    palette='colorblind'
)
axes[0].set_title('Donor Site (Upstream Exon)', fontsize=14, fontweight='bold')
axes[0].set_xlabel('Position relative to splice site', fontsize=12)
axes[0].set_ylabel('Mean Architecture Fraction', fontsize=12)
axes[0].set_xticks(range(-20, 0, 2))
axes[0].axvline(x=-0.5, color='red', linestyle='--', alpha=0.7, label='Splice Site')

# Plot 2: Acceptor Site (+1 to +20)
acceptor_df = df[df['site'] == 'acceptor'].copy()
acceptor_df = acceptor_df.sort_values('position')

sns.lineplot(
    data=acceptor_df,
    x='position',
    y='mean_arch_frac',
    hue='dominant_codon',
    hue_order=codon_order,
    marker='o',
    ax=axes[1],
    palette='colorblind'
)
axes[1].set_title('Acceptor Site (Downstream Exon)', fontsize=14, fontweight='bold')
axes[1].set_xlabel('Position relative to splice site', fontsize=12)
axes[1].set_ylabel('')
axes[1].set_xticks(range(1, 21, 2))
axes[1].axvline(x=0.5, color='red', linestyle='--', alpha=0.7)

# Global formatting
plt.suptitle('Splice Site Architecture Fraction by Codon Position', fontsize=16, fontweight='bold', y=1.05)
plt.tight_layout()

# Save high-resolution figure
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"Plot saved to {OUTPUT_PLOT}")
