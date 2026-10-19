import pandas as pd

INPUT_FILE = "data/processed/correct_phyloP_profiles.tsv"

df = pd.read_csv(INPUT_FILE, sep="\t")

row = df[
    (df["gene_name"] == "DNMBP") &
    (df["strand"] == "-")
].iloc[0]

donor = [
    float(x) if x != "nan" else None
    for x in row["donor_profile"].split(",")
]

acceptor = [
    float(x) if x != "nan" else None
    for x in row["acceptor_profile"].split(",")
]

print("========================================")
print("CORRECT PHYLOP PROFILE CHECK")
print("========================================")

print("Gene:", row["gene_name"])
print("Junction:", row["junction_id"])
print("Strand:", row["strand"])

print("\nDONOR PROFILE")
print("----------------------------------------")

positions = list(range(-20, 0)) + list(range(1, 21))

for position, score in zip(positions, donor):
    print(f"{position:>3}: {score:.4f}")

print("\nExpected donor motif:")
print("+1 = G")
print("+2 = T")

print("\nDONOR MOTIF PHYLOP:")
print("+1:", donor[20])
print("+2:", donor[21])


print("\nACCEPTOR PROFILE")
print("----------------------------------------")

for position, score in zip(positions, acceptor):
    print(f"{position:>3}: {score:.4f}")

print("\nExpected acceptor motif:")
print("-2 = A")
print("-1 = G")

print("\nACCEPTOR MOTIF PHYLOP:")
print("-2:", acceptor[18])
print("-1:", acceptor[19])
