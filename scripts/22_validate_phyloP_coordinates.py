import pandas as pd
import pysam
import pyBigWig

INPUT_FILE = "data/processed/splice_site_features.tsv"
GENOME = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"
PHYLOP = "data/reference/conservation/hg38.phyloP100way.bw"

df = pd.read_csv(INPUT_FILE, sep="\t")

# Select the same type of negative-strand canonical junction
row = df[
    (df["gene_name"] == "DNMBP") &
    (df["strand"] == "-") &
    (df["donor_motif"] == "GT") &
    (df["acceptor_motif"] == "AG")
].iloc[0]

chrom = row["chromosome"]
strand = row["strand"]

intron_start = int(row["intron_start"])
intron_end = int(row["intron_end"])

genome = pysam.FastaFile(GENOME)
bw = pyBigWig.open(PHYLOP)

print("========================================")
print("PHYLOP COORDINATE VALIDATION")
print("========================================")

print("Gene:", row["gene_name"])
print("Chromosome:", chrom)
print("Strand:", strand)

print("Intron:", intron_start, "-", intron_end)

print("\nExpected donor:", row["donor_motif"])
print("Expected acceptor:", row["acceptor_motif"])


# ============================================================
# DONOR
# ============================================================

print("\n========================================")
print("DONOR")
print("========================================")

# For negative strand:
# donor GT corresponds to genomic AC at:
# intron_end-1 and intron_end

donor_genomic = genome.fetch(
    chrom,
    intron_end - 2,
    intron_end
).upper()

print("Genomic donor sequence:", donor_genomic)

complement = str.maketrans(
    "ACGT",
    "TGCA"
)

donor_transcript = donor_genomic.translate(
    complement
)[::-1]

print("Transcript donor sequence:", donor_transcript)

print("\nGenomic donor positions:")

for pos in range(intron_end - 5, intron_end + 6):

    base = genome.fetch(
        chrom,
        pos - 1,
        pos
    ).upper()

    score = bw.values(
        chrom,
        pos - 1,
        pos
    )[0]

    print(
        f"{pos}: {base}   phyloP={score:.4f}"
    )


# ============================================================
# ACCEPTOR
# ============================================================

print("\n========================================")
print("ACCEPTOR")
print("========================================")

acceptor_genomic = genome.fetch(
    chrom,
    intron_start - 1,
    intron_start + 1
).upper()

print("Genomic acceptor sequence:", acceptor_genomic)

acceptor_transcript = acceptor_genomic.translate(
    complement
)[::-1]

print(
    "Transcript acceptor sequence:",
    acceptor_transcript
)

print("\nGenomic acceptor positions:")

for pos in range(intron_start - 5, intron_start + 6):

    base = genome.fetch(
        chrom,
        pos - 1,
        pos
    ).upper()

    score = bw.values(
        chrom,
        pos - 1,
        pos
    )[0]

    print(
        f"{pos}: {base}   phyloP={score:.4f}"
    )


genome.close()
