import pysam

genome = "data/reference/GRCh38.primary_assembly.genome.fa.bgz"

fasta = pysam.FastaFile(genome)

chrom = "chrX"

positions = range(100633395, 100633416)

print("Position    Base")
print("----------------")

for pos in positions:
    base = fasta.fetch(
        chrom,
        pos - 1,
        pos
    )
    print(f"{pos}    {base}")

fasta.close()
