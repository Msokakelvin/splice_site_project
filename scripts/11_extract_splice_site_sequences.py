import pandas as pd
import pysam

junction_file = "data/processed/alternative_junction_summary.tsv"
genome = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"
output_file = "data/processed/splice_site_sequences.tsv"

WINDOW = 20


def reverse_complement(seq):
    complement = str.maketrans(
        "ACGTN",
        "TGCAN"
    )
    return seq.translate(complement)[::-1]


print("Opening reference genome...")

fasta = pysam.FastaFile(genome)

print("Reading junction table...")

df = pd.read_csv(
    junction_file,
    sep="\t"
)

print("Junctions:", len(df))

results = []

for i, row in df.iterrows():

    chromosome = row["chromosome"]
    strand = row["strand"]

    donor = int(row["donor"])
    acceptor = int(row["acceptor"])

    # ==================================================
    # DONOR SITE
    # ==================================================

    if strand == "+":

        # Donor is exon end.
        # First two intronic bases are immediately
        # after the exon.
        donor_motif_start = donor + 1

        donor_seq = fasta.fetch(
            chromosome,
            donor - WINDOW,
            donor + WINDOW + 2
        )

        donor_motif = fasta.fetch(
            chromosome,
            donor,
            donor + 2
        ).upper()

    else:

        # Donor is exon start.
        # First two intronic bases are immediately
        # before the exon.
        donor_motif_start = donor - 2

        donor_seq = fasta.fetch(
            chromosome,
            donor - WINDOW - 2,
            donor + WINDOW
        )

        donor_motif = fasta.fetch(
            chromosome,
            donor - 2,
            donor
        ).upper()

        donor_motif = reverse_complement(donor_motif)

        donor_seq = reverse_complement(
            donor_seq
        )

    # ==================================================
    # ACCEPTOR SITE
    # ==================================================

    if strand == "+":

        # Acceptor is exon start.
        # Last two intronic bases are immediately
        # before the exon.
        acceptor_seq = fasta.fetch(
            chromosome,
            acceptor - WINDOW - 2,
            acceptor + WINDOW
        )

        acceptor_motif = fasta.fetch(
            chromosome,
            acceptor - 2,
            acceptor
        ).upper()

    else:

        # Acceptor is exon end.
        # Last two intronic bases are immediately
        # after the exon in genomic coordinates.
        acceptor_seq = fasta.fetch(
            chromosome,
            acceptor - WINDOW,
            acceptor + WINDOW + 2
        )

        acceptor_motif = fasta.fetch(
            chromosome,
            acceptor,
            acceptor + 2
        ).upper()

        acceptor_motif = reverse_complement(
            acceptor_motif
        )

        acceptor_seq = reverse_complement(
            acceptor_seq
        )

    # ==================================================
    # CLASSIFY SPLICE SITE
    # ==================================================

    if donor_motif == "GT" and acceptor_motif == "AG":

        splice_class = "canonical_GT_AG"

    elif donor_motif == "GC" and acceptor_motif == "AG":

        splice_class = "noncanonical_GC_AG"

    elif donor_motif == "AT" and acceptor_motif == "AC":

        splice_class = "noncanonical_AT_AC"

    else:

        splice_class = "other"

    results.append([
        row["junction_id"],
        row["gene_id"],
        chromosome,
        strand,
        donor,
        acceptor,
        row["events"],
        row["transcripts"],
        row["roles"],
        row["role_class"],
        donor_seq.upper(),
        acceptor_seq.upper(),
        donor_motif,
        acceptor_motif,
        splice_class
    ])

    if (i + 1) % 10000 == 0:

        print(
            f"Processed {i + 1:,}/{len(df):,}"
        )


result = pd.DataFrame(
    results,
    columns=[
        "junction_id",
        "gene_id",
        "chromosome",
        "strand",
        "donor",
        "acceptor",
        "events",
        "transcripts",
        "roles",
        "role_class",
        "donor_sequence",
        "acceptor_sequence",
        "donor_motif",
        "acceptor_motif",
        "splice_class"
    ]
)

result.to_csv(
    output_file,
    sep="\t",
    index=False
)

fasta.close()

print()
print("Finished.")
print("Sequences:", len(result))

print()
print("Splice classes:")
print(
    result["splice_class"].value_counts()
)

print()
print("Output:", output_file)
