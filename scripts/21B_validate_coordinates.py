import pandas as pd
import pysam

INPUT_FILE = "data/processed/splice_site_features.tsv"
GENOME = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"

df = pd.read_csv(INPUT_FILE, sep="\t")
genome = pysam.FastaFile(GENOME)

complement = str.maketrans(
    "ACGT",
    "TGCA"
)

checked = 0
correct = 0

for _, row in df.head(10000).iterrows():

    chrom = row["chromosome"]
    strand = row["strand"]

    intron_start = int(row["intron_start"])
    intron_end = int(row["intron_end"])

    expected_donor = row["donor_motif"]
    expected_acceptor = row["acceptor_motif"]

    # ==========================================
    # POSITIVE STRAND
    # ==========================================

    if strand == "+":

        donor_seq = genome.fetch(
            chrom,
            intron_start - 1,
            intron_start + 1
        ).upper()

        acceptor_seq = genome.fetch(
            chrom,
            intron_end - 2,
            intron_end
        ).upper()

    # ==========================================
    # NEGATIVE STRAND
    # ==========================================

    else:

        donor_genomic = genome.fetch(
            chrom,
            intron_end - 2,
            intron_end
        ).upper()

        acceptor_genomic = genome.fetch(
            chrom,
            intron_start - 1,
            intron_start + 1
        ).upper()

        donor_seq = donor_genomic.translate(
            complement
        )[::-1]

        acceptor_seq = acceptor_genomic.translate(
            complement
        )[::-1]

    donor_ok = donor_seq == expected_donor
    acceptor_ok = acceptor_seq == expected_acceptor

    if donor_ok and acceptor_ok:
        correct += 1

    checked += 1

    # Show only errors
    if not (donor_ok and acceptor_ok) and checked <= 20:

        print("\nERROR")
        print("Gene:", row["gene_name"])
        print("Strand:", strand)
        print("Intron:", intron_start, "-", intron_end)

        print("Expected donor:", expected_donor)
        print("Observed donor:", donor_seq)

        print("Expected acceptor:", expected_acceptor)
        print("Observed acceptor:", acceptor_seq)


genome.close()

print("\n========================================")
print("DEFINITIVE MOTIF VALIDATION")
print("========================================")

print("Checked:", checked)
print("Correct:", correct)

print(
    "Validation rate:",
    f"{100 * correct / checked:.2f}%"
)
