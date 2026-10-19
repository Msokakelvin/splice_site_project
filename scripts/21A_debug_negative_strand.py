import pandas as pd
import pysam

INPUT_FILE = "data/processed/splice_site_features.tsv"
GENOME = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"

df = pd.read_csv(INPUT_FILE, sep="\t")
genome = pysam.FastaFile(GENOME)

# Select the first negative-strand canonical GT-AG junction
row = df[
    (df["strand"] == "-") &
    (df["donor_motif"] == "GT") &
    (df["acceptor_motif"] == "AG")
].iloc[0]

chrom = row["chromosome"]
intron_start = int(row["intron_start"])
intron_end = int(row["intron_end"])

print("========================================")
print("NEGATIVE-STRAND SPLICE SITE DEBUG")
print("========================================")

print("Gene:", row["gene_name"])
print("Gene ID:", row["gene_id"])
print("Chromosome:", chrom)
print("Strand:", row["strand"])

print("\nIntron coordinates:")
print("intron_start:", intron_start)
print("intron_end:", intron_end)

print("\nExpected:")
print("donor:", row["donor_motif"])
print("acceptor:", row["acceptor_motif"])

# ------------------------------------------------
# Show sequence around the HIGH-coordinate end
# ------------------------------------------------

print("\n----------------------------------------")
print("HIGH-COORDINATE END OF INTRON")
print("----------------------------------------")

start = intron_end - 10
end = intron_end + 10

seq = genome.fetch(
    chrom,
    start - 1,
    end
).upper()

print("Genomic positions:", start, "to", end)
print("Sequence:", seq)

# The last two bases inside the intron
donor_genomic = genome.fetch(
    chrom,
    intron_end - 2,
    intron_end
).upper()

print("\nGenomic bases at intron_end-1 to intron_end:")
print(donor_genomic)

complement = str.maketrans(
    "ACGT",
    "TGCA"
)

donor_transcript = donor_genomic.translate(
    complement
)[::-1]

print("Reverse complement:")
print(donor_transcript)

# ------------------------------------------------
# Try the neighboring possibilities
# ------------------------------------------------

print("\n========================================")
print("TESTING NEIGHBORING TWO-BASE WINDOWS")
print("========================================")

for offset in range(-5, 6):

    pos1 = intron_end + offset
    pos2 = pos1 + 1

    genomic = genome.fetch(
        chrom,
        pos1 - 1,
        pos2
    ).upper()

    transcript = genomic.translate(
        complement
    )[::-1]

    print(
        f"Genomic positions {pos1}-{pos2}: "
        f"{genomic} -> transcript orientation {transcript}"
    )

# ------------------------------------------------
# Acceptor check
# ------------------------------------------------

print("\n========================================")
print("ACCEPTOR CHECK")
print("========================================")

acceptor_genomic = genome.fetch(
    chrom,
    intron_start - 1,
    intron_start + 1
).upper()

acceptor_transcript = acceptor_genomic.translate(
    complement
)[::-1]

print("Genomic sequence:", acceptor_genomic)
print("Transcript orientation:", acceptor_transcript)
print("Expected:", row["acceptor_motif"])

genome.close()
