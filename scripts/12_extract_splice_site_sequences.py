import pandas as pd
import pysam

intron_file = "data/processed/intron_coordinates.tsv"
output_file = "data/processed/splice_site_sequences.tsv"

genome = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"

WINDOW = 20


def reverse_complement(seq):

    complement = str.maketrans(
        "ACGTN",
        "TGCAN"
    )

    return seq.translate(complement)[::-1]


print("Opening reference genome...")

fasta = pysam.FastaFile(genome)


print("Reading intron table...")

df = pd.read_csv(
    intron_file,
    sep="\t"
)

print("Introns:", len(df))


results = []


for i, row in df.iterrows():

    chromosome = row["chromosome"]
    strand = row["strand"]

    intron_start = int(row["intron_start"])
    intron_end = int(row["intron_end"])


    # ==================================================
    # DONOR AND ACCEPTOR MOTIFS
    # ==================================================

    if strand == "+":

        # 5' end of intron = donor
        donor_motif = fasta.fetch(
            chromosome,
            intron_start - 1,
            intron_start + 1
        ).upper()

        # 3' end of intron = acceptor
        acceptor_motif = fasta.fetch(
            chromosome,
            intron_end - 2,
            intron_end
        ).upper()


        # Sequence windows

        donor_sequence = fasta.fetch(
            chromosome,
            intron_start - WINDOW - 1,
            intron_start + WINDOW + 1
        ).upper()

        acceptor_sequence = fasta.fetch(
            chromosome,
            intron_end - WINDOW,
            intron_end + WINDOW
        ).upper()


    else:

        # On minus strand:
        # genomic intron end = transcript 5' donor
        donor_motif = fasta.fetch(
            chromosome,
            intron_end - 2,
            intron_end
        ).upper()

        # genomic intron start = transcript 3' acceptor
        acceptor_motif = fasta.fetch(
            chromosome,
            intron_start - 1,
            intron_start + 1
        ).upper()


        # Convert to transcript orientation

        donor_motif = reverse_complement(
            donor_motif
        )

        acceptor_motif = reverse_complement(
            acceptor_motif
        )


        donor_sequence = fasta.fetch(
            chromosome,
            intron_end - WINDOW - 1,
            intron_end + WINDOW
        ).upper()

        acceptor_sequence = fasta.fetch(
            chromosome,
            intron_start - WINDOW - 1,
            intron_start + WINDOW
        ).upper()


        donor_sequence = reverse_complement(
            donor_sequence
        )

        acceptor_sequence = reverse_complement(
            acceptor_sequence
        )


    # ==================================================
    # CLASSIFY SPLICE TYPE
    # ==================================================

    if donor_motif == "GT" and acceptor_motif == "AG":

        splice_class = "canonical_GT_AG"

    elif donor_motif == "GC" and acceptor_motif == "AG":

        splice_class = "noncanonical_GC_AG"

    elif donor_motif == "AT" and acceptor_motif == "AC":

        splice_class = "noncanonical_AT_AC"

    else:

        splice_class = "other"


    # ==================================================
    # CREATE UNIQUE INTRON ID
    # ==================================================

    intron_id = (
        chromosome
        + ":"
        + str(intron_start)
        + "-"
        + str(intron_end)
        + ":"
        + strand
    )


    results.append([

        intron_id,

        row["gene_id"],

        row["gene_name"],

        row["transcript_id"],

        chromosome,

        strand,

        intron_start,

        intron_end,

        donor_motif,

        acceptor_motif,

        splice_class,

        donor_sequence,

        acceptor_sequence

    ])


    if (i + 1) % 10000 == 0:

        print(
            f"Processed "
            f"{i + 1:,}/"
            f"{len(df):,}"
        )


# ==================================================
# CREATE OUTPUT TABLE
# ==================================================

result = pd.DataFrame(

    results,

    columns=[

        "intron_id",

        "gene_id",

        "gene_name",

        "transcript_id",

        "chromosome",

        "strand",

        "intron_start",

        "intron_end",

        "donor_motif",

        "acceptor_motif",

        "splice_class",

        "donor_sequence",

        "acceptor_sequence"

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
print("Intron sequences:", len(result))

print()
print("Splice classes:")

print(
    result["splice_class"].value_counts()
)

print()
print("Output:", output_file)
