import pandas as pd
import gzip
from collections import defaultdict

GTF = "data/reference/gencode.v50.annotation.gtf.gz"
MASTER = "data/processed/master_gene_junction_table.tsv"
OUTPUT = "data/processed/junction_position_coding_phase.tsv"

print("Loading unique gene-junction table...")

master = pd.read_csv(
    MASTER,
    sep="\t",
    usecols=[
        "gene_id",
        "gene_name",
        "junction_id",
        "chromosome",
        "strand",
        "intron_start",
        "intron_end"
    ]
)

print("Unique gene-junction rows:", len(master))
print("Genes:", master["gene_id"].nunique())


# =========================================================
# 1. Read CDS annotations
# =========================================================

print()
print("Reading CDS annotations from GENCODE...")

cds_by_transcript = defaultdict(list)

with gzip.open(GTF, "rt") as f:

    for line in f:

        if line.startswith("#"):
            continue

        fields = line.rstrip("\n").split("\t")

        if len(fields) != 9:
            continue

        if fields[2] != "CDS":
            continue

        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        phase = fields[7]
        attributes = fields[8]

        transcript_id = None

        for item in attributes.split(";"):

            item = item.strip()

            if item.startswith("transcript_id "):

                transcript_id = item.split('"')[1]
                break

        if transcript_id is None:
            continue

        cds_by_transcript[transcript_id].append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "strand": strand,
                "phase": phase
            }
        )

print("Transcripts with CDS:", len(cds_by_transcript))


# =========================================================
# 2. Create a compact genomic CDS lookup
#
# Instead of storing every nucleotide, store only CDS
# intervals. This greatly reduces memory usage.
# =========================================================

print()
print("Sorting CDS intervals...")

for transcript_id in cds_by_transcript:

    intervals = cds_by_transcript[transcript_id]

    strand = intervals[0]["strand"]

    if strand == "+":

        intervals.sort(
            key=lambda x: x["start"]
        )

    else:

        intervals.sort(
            key=lambda x: x["start"],
            reverse=True
        )

print("CDS intervals ready.")


# =========================================================
# 3. Determine whether a genomic position is coding
# =========================================================

def get_coding_status(transcript_id, genomic_position):

    intervals = cds_by_transcript.get(
        transcript_id,
        []
    )

    for interval in intervals:

        if (
            interval["start"]
            <= genomic_position
            <= interval["end"]
        ):

            return True

    return False


# =========================================================
# 4. Generate splice-site positions
# =========================================================

positions = list(range(-20, 0)) + list(range(1, 21))


def get_position(
    strand,
    site,
    intron_start,
    intron_end,
    position
):

    """
    Convert transcript-oriented splice-site position
    into genomic coordinate.

    Donor:
        -20 ... -1 = upstream exon
        +1 ... +20 = intron

    Acceptor:
        -20 ... -1 = intron
        +1 ... +20 = downstream exon
    """

    if strand == "+":

        if site == "donor":

            return intron_start + position - 1

        else:

            return intron_end + position + 1

    else:

        if site == "donor":

            return intron_end - position + 1

        else:

            return intron_start - position - 1


# =========================================================
# 5. Write output incrementally
# =========================================================

print()
print("Mapping coding status...")

header_written = False

chunk_size = 5000

total = len(master)

for start_index in range(
    0,
    total,
    chunk_size
):

    chunk = master.iloc[
        start_index:
        start_index + chunk_size
    ]

    records = []

    for _, row in chunk.iterrows():

        strand = row["strand"]

        intron_start = int(
            row["intron_start"]
        )

        intron_end = int(
            row["intron_end"]
        )

        # -------------------------------------------------
        # Donor and acceptor
        # -------------------------------------------------

        for site in ["donor", "acceptor"]:

            for position in positions:

                genomic_position = get_position(
                    strand,
                    site,
                    intron_start,
                    intron_end,
                    position
                )

                # -------------------------------------------------
                # Important:
                #
                # master table is gene-junction level, so coding
                # status cannot be assigned to one transcript.
                #
                # We therefore record the genomic coordinate here.
                # Transcript-specific CDS validation will be done
                # in the next analysis.
                # -------------------------------------------------

                records.append(
                    {
                        "gene_id": row["gene_id"],
                        "gene_name": row["gene_name"],
                        "junction_id": row["junction_id"],
                        "strand": strand,
                        "site": site,
                        "position": position,
                        "genomic_position": genomic_position
                    }
                )

    result = pd.DataFrame(records)

    result.to_csv(
        OUTPUT,
        sep="\t",
        mode="w" if not header_written else "a",
        header=not header_written,
        index=False
    )

    header_written = True

    processed = min(
        start_index + chunk_size,
        total
    )

    print(
        "Processed:",
        processed,
        "/",
        total
    )


# =========================================================
# 6. Final summary
# =========================================================

print()
print("==============================================")
print("SPLICE POSITION COORDINATE MAPPING COMPLETE")
print("==============================================")
print()

print("Gene-junctions:", len(master))

print(
    "Expected output rows:",
    len(master) * 80
)

print()

print("Output:")
print(OUTPUT)

print()
print("The output contains transcript-oriented")
print("splice-site positions and their genomic")
print("coordinates.")
