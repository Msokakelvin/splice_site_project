import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

INPUT_FILE = "data/processed/final_junction_analysis_table.tsv"
OUTPUT_PLOT = "results/constitutive_vs_alternative_conservation.png"
OUTPUT_STATS = "data/processed/constitutive_vs_alternative_stats.tsv"

print("=== Script 42: Constitutive vs. Alternative Conservation ===")

# 1. Load data
print("Loading final table...")
df = pd.read_csv(INPUT_FILE, sep="\t")

# 2. Define Constitutive vs. Alternative
# We assume 'constitutive' is explicitly labeled in splice_class. 
# If not, we can define it as architecture_fraction == 1.0. Let's use splice_class first.
if 'constitutive' in df['splice_class'].str.lower().values:
    df['junction_type'] = df['splice_class'].str.lower().apply(
        lambda x: 'Constitutive' if 'constitutive' in x else 'Alternative'
    )
else:
    # Fallback: define constitutive as architecture_fraction == 1.0 (or > 0.95)
    df['junction_type'] = df['architecture_fraction'].apply(
        lambda x: 'Constitutive' if x >= 0.95 else 'Alternative'
    )

print(f"Found {df['junction_type'].value_counts().to_dict()} junction types.")

# 3. Statistical Test (Mann-Whitney U test, since architecture fraction is not normally distributed)
constitutive_scores = df[df['junction_type'] == 'Constitutive']['architecture_fraction']
alternative_scores = df[df['junction_type'] == 'Alternative']['architecture_fraction']

stat, p_value = stats.mannwhitneyu(constitutive_scores, alternative_scores, alternative='greater')

print("\n--- Statistical Test (Mann-Whitney U) ---")
print(f"Constitutive mean architecture fraction: {constitutive_scores.mean():.4f}")
print(f"Alternative mean architecture fraction:  {alternative_scores.mean():.4f}")
print(f"P-value: {p_value:.2e}")
if p_value < 0.001:
    print("Result: Constitutive junctions are significantly more conserved than alternative junctions. ✅")

# Save stats
stats_df = pd.DataFrame({
    'junction_type': ['Constitutive', 'Alternative'],
    'mean_arch_frac': [constitutive_scores.mean(), alternative_scores.mean()],
    'median_arch_frac': [constitutive_scores.median(), alternative_scores.median()],
    'count': [len(constitutive_scores), len(alternative_scores)]
})
stats_df.to_csv(OUTPUT_STATS, sep="\t", index=False)
print(f"\nStats saved to {OUTPUT_STATS}")

# 4. Plotting
print("\nGenerating plot...")
sns.set_theme(style="whitegrid", font_scale=1.2)
plt.figure(figsize=(8, 6))

# Use a violin plot with a boxplot inside to show distribution and median
sns.violinplot(
    data=df,
    x='junction_type',
    y='architecture_fraction',
    palette='Set2',
    inner='box',
    cut=0
)

plt.title('Sequence Conservation: Constitutive vs. Alternative Splice Sites', fontsize=14, fontweight='bold')
plt.xlabel('Junction Type', fontsize=12)
plt.ylabel('Architecture Fraction', fontsize=12)

# Add significance star
max_y = df['architecture_fraction'].max()
plt.text(0.5, max_y * 0.95, f'p < 0.001' if p_value < 0.001 else f'p = {p_value:.2e}', 
         ha='center', va='center', fontsize=14, fontweight='bold', color='red')

plt.tight_layout()
plt.savefig(OUTPUT_PLOT, dpi=300, bbox_inches='tight')
print(f"Plot saved to {OUTPUT_PLOT}")
