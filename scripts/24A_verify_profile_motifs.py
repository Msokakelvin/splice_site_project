import pandas as pd
import pysam

FEATURE_FILE = "data/processed/splice_site_features.tsv"
PROFILE_FILE = "data/processed/correct_phyloP_profiles.tsv"
GENOME_FILE = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"

TARGET_INTRON = "chr10:99877337-99879810:-"


def reverse_complement(seq):
    complement = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan"
    )
    return seq.translate(complement)[::-1].upper()


# --------------------------------------------------
# Load files
# --------------------------------------------------

features = pd.read_csv(
    FEATURE_FILE,
    sep="\t"
)

profiles = pd.read_csv(
    PROFILE_FILE,
    sep="\t"
)

genome = pysam.FastaFile(
    GENOME_FILE
)


# --------------------------------------------------
# Find target intron in feature table
# --------------------------------------------------

matches = features[
    features["intron_id"] == TARGET_INTRON
]

if matches.empty:
    raise ValueError(
        f"Target intron not found: {TARGET_INTRON}"
    )

feature_row = matches.iloc[0]


# --------------------------------------------------
# Extract information
# --------------------------------------------------

gene = feature_row["gene_name"]
chrom = feature_row["chromosome"]
strand = feature_row["strand"]

intron_start = int(
    feature_row["intron_start"]
)

intron_end = int(
    feature_row["intron_end"]
)


print()
print("PROFILE MOTIF VERIFICATION")
print("=" * 60)

print(f"Gene:       {gene}")
print(f"Intron:     {TARGET_INTRON}")
print(f"Chromosome: {chrom}")
print(f"Strand:     {strand}")

print(
    f"Intron:     {intron_start}-{intron_end}"
)

print()
print("Expected:")
print("Donor:      GT")
print("Acceptor:   AG")


# --------------------------------------------------
# Extract genomic splice motifs
# --------------------------------------------------

if strand == "+":

    donor_genomic = genome.fetch(
        chrom,
        intron_start - 1,
        intron_start + 1
    ).upper()

    acceptor_genomic = genome.fetch(
        chrom,
        intron_end - 2,
        intron_end
    ).upper()

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


# --------------------------------------------------
# Convert to transcript orientation
# --------------------------------------------------

if strand == "-":

    donor_transcript = reverse_complement(
        donor_genomic
    )

    acceptor_transcript = reverse_complement(
        acceptor_genomic
    )

else:

    donor_transcript = donor_genomic
    acceptor_transcript = acceptor_genomic


# --------------------------------------------------
# Display motif verification
# --------------------------------------------------

print()
print("DONOR")
print("-" * 60)

print(
    f"Genomic motif:    {donor_genomic}"
)

print(
    f"Transcript motif: {donor_transcript}"
)


print()
print("ACCEPTOR")
print("-" * 60)

print(
    f"Genomic motif:    {acceptor_genomic}"
)

print(
    f"Transcript motif: {acceptor_transcript}"
)


# --------------------------------------------------
# Find corresponding profile using coordinates
# --------------------------------------------------

profile_matches = profiles[
    (profiles["intron_start"] == intron_start)
    &
    (profiles["intron_end"] == intron_end)
    &
    (profiles["strand"] == strand)
]


if profile_matches.empty:

    print()
    print(
        "WARNING: Matching phyloP profile not found."
    )

else:

    profile_row = profile_matches.iloc[0]

    print()
    print("PHYLOP PROFILE")
    print("-" * 60)

    print(
        f"Profile junction ID: "
        f"{profile_row['junction_id']}"
    )

    print()

    # Split the profile strings into numbers
    donor_values = [
        float(x)
        for x in str(
            profile_row["donor_profile"]
        ).split(",")
    ]

    acceptor_values = [
        float(x)
        for x in str(
            profile_row["acceptor_profile"]
        ).split(",")
    ]


    print(
        f"Donor profile length: "
        f"{len(donor_values)}"
    )

    print(
        f"Acceptor profile length: "
        f"{len(acceptor_values)}"
    )

    print()

    # Profile positions:
    #
    # Donor:
    # -20 ... -1 | +1 ... +20
    #
    # Acceptor:
    # -20 ... -2 -1 | +1 ... +20

    donor_plus1 = donor_values[20]
    donor_plus2 = donor_values[21]

    acceptor_minus2 = acceptor_values[18]
    acceptor_minus1 = acceptor_values[19]


    print("PHYLOP VALUES AT SPLICE MOTIF")
    print("-" * 60)

    print(
        f"Donor +1:     {donor_plus1:.3f}"
    )

    print(
        f"Donor +2:     {donor_plus2:.3f}"
    )

    print(
        f"Acceptor -2:  {acceptor_minus2:.3f}"
    )

    print(
        f"Acceptor -1:  {acceptor_minus1:.3f}"
    )


genome.close()


print()
print("=" * 60)
print("Verification complete.")
