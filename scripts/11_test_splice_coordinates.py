import pysam

genome = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"

fasta = pysam.FastaFile(genome)


def reverse_complement(seq):
    return seq.translate(
        str.maketrans("ACGTN", "TGCAN")
    )[::-1]


chrom = "chrX"
strand = "-"
donor = 100633405
acceptor = 100632568


print("Testing known splice junction")
print("--------------------------------")
print("Chromosome:", chrom)
print("Strand:", strand)
print("Donor boundary:", donor)
print("Acceptor boundary:", acceptor)


# -------------------------
# DONOR
# -------------------------

print("\nDONOR")

if strand == "+":

    # donor is exon end
    genomic_donor = fasta.fetch(
        chrom,
        donor,
        donor + 2
    )

else:

    # donor is exon start
    # two genomic bases immediately before exon
    genomic_donor = fasta.fetch(
        chrom,
        donor - 2,
        donor
    )

print("Genomic donor bases:", genomic_donor)

if strand == "-":
    donor_motif = reverse_complement(genomic_donor)
else:
    donor_motif = genomic_donor

print("Transcript donor motif:", donor_motif)


# -------------------------
# ACCEPTOR
# -------------------------

print("\nACCEPTOR")

if strand == "+":

    # acceptor is exon start
    # two genomic bases immediately before exon
    genomic_acceptor = fasta.fetch(
        chrom,
        acceptor - 2,
        acceptor
    )

else:

    # acceptor is exon end
    genomic_acceptor = fasta.fetch(
        chrom,
        acceptor,
        acceptor + 2
    )

print("Genomic acceptor bases:", genomic_acceptor)

if strand == "-":
    acceptor_motif = reverse_complement(genomic_acceptor)
else:
    acceptor_motif = genomic_acceptor

print("Transcript acceptor motif:", acceptor_motif)


# -------------------------
# CLASSIFICATION
# -------------------------

print("\nRESULT")

if donor_motif == "GT" and acceptor_motif == "AG":
    print("Canonical splice site: GT-AG")

elif donor_motif == "GC" and acceptor_motif == "AG":
    print("Noncanonical splice site: GC-AG")

elif donor_motif == "AT" and acceptor_motif == "AC":
    print("Noncanonical splice site: AT-AC")

else:
    print("Other:", donor_motif, "-", acceptor_motif)


fasta.close()
