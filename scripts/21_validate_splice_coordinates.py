import pandas as pd
import pysam

INPUT_FILE = "data/processed/splice_site_features.tsv"
GENOME = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"

df = pd.read_csv(INPUT_FILE, sep="\t")
genome = pysam.FastaFile(GENOME)

print("Checking splice motifs...")

checked = 0
correct = 0

for _, row in df.head(10000).iterrows():

    chrom = row["chromosome"]
    strand = row["strand"]

    intron_start = int(row["intron_start"])
    intron_end = int(row["intron_end"])

    if strand == "+":

        # Donor = first 2 bases of intron
        donor_genomic = genome.fetch(
    chrom,
    intron_end - 2,
    intron_end
    ).upper()

        # Acceptor = last 2 bases of intron
        acceptor_seq = genome.fetch(
            chrom,
            intron_end - 2,
            intron_end
        ).upper()

    else:

        # Extract genomic sequences first
        donor_genomic = genome.fetch(
            chrom,
            intron_end - 1,
            intron_end + 1
        ).upper()

        acceptor_genomic = genome.fetch(
            chrom,
            intron_start - 1,
            intron_start + 1
        ).upper()

        complement = str.maketrans(
            "ACGT",
            "TGCA"
        )

        # Convert genomic orientation to transcript orientation
        donor_seq = donor_genomic.translate(
            complement
        )[::-1]

        acceptor_seq = acceptor_genomic.translate(
            complement
        )[::-1]

    expected_donor = row["donor_motif"]
    expected_acceptor = row["acceptor_motif"]

    donor_ok = donor_seq == expected_donor
    acceptor_ok = acceptor_seq == expected_acceptor

    if donor_ok and acceptor_ok:
        correct += 1

    checked += 1

    if checked <= 20:

        print("\nGene:", row["gene_name"])
        print("Strand:", strand)

        print("Expected donor:", expected_donor)
        print("Observed donor:", donor_seq)

        print("Expected acceptor:", expected_acceptor)
        print("Observed acceptor:", acceptor_seq)

        print("Correct:", donor_ok and acceptor_ok)


genome.close()


print("\n========================================")
print("MOTIF COORDINATE VALIDATION")
print("========================================")

print("Checked:", checked)
print("Correct donor + acceptor:", correct)

print(
    "Validation rate:",
    f"{100 * correct / checked:.2f}%"
)
