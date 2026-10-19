#!/usr/bin/env python3

import gzip
import bisect
from collections import defaultdict

import numpy as np
import pandas as pd


# ============================================================
# FILES
# ============================================================

GTF = "data/reference/gencode.v50.annotation.gtf.gz"

SPLICE_JUNCTIONS = "data/processed/splice_junctions.tsv"

MASTER = "data/processed/master_gene_junction_table.tsv"

OUTPUT = "data/processed/junction_coding_context.tsv"


# ============================================================
# SETTINGS
# ============================================================

CHUNK_SIZE = 100_000


# ============================================================
# FUNCTION 1
# LOAD CDS ANNOTATION FROM GENCODE
# ============================================================

def load_cds_annotation():

    print()
    print("=" * 70)
    print("LOADING CDS ANNOTATION FROM GENCODE")
    print("=" * 70)

    cds_by_transcript = defaultdict(list)

    phase_counts = {
        "0": 0,
        "1": 0,
        "2": 0,
        "unknown": 0
    }

    with gzip.open(GTF, "rt") as handle:

        for line in handle:

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            if len(fields) != 9:
                continue

            feature = fields[2]

            if feature != "CDS":
                continue

            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            phase = fields[7]
            attributes = fields[8]

            # ------------------------------------------------
            # Extract transcript_id
            # ------------------------------------------------

            transcript_id = None

            for item in attributes.split(";"):

                item = item.strip()

                if item.startswith("transcript_id "):

                    parts = item.split('"')

                    if len(parts) >= 2:
                        transcript_id = parts[1]

                    break

            if transcript_id is None:
                continue

            # ------------------------------------------------
            # CDS phase
            # ------------------------------------------------

            if phase in ("0", "1", "2"):

                phase_int = int(phase)

                phase_counts[phase] += 1

            else:

                phase_int = None

                phase_counts["unknown"] += 1

            # ------------------------------------------------
            # Store CDS interval
            # ------------------------------------------------

            cds_by_transcript[transcript_id].append(
                (
                    start,
                    end,
                    strand,
                    phase_int
                )
            )

    # --------------------------------------------------------
    # Sort CDS intervals by genomic start
    # --------------------------------------------------------

    for transcript_id in cds_by_transcript:

        cds_by_transcript[transcript_id].sort(
            key=lambda x: x[0]
        )

    # --------------------------------------------------------
    # Store genomic starts separately for fast searching
    # --------------------------------------------------------

    cds_starts = {}

    for transcript_id, intervals in cds_by_transcript.items():

        cds_starts[transcript_id] = [
            interval[0]
            for interval in intervals
        ]

    print(
        f"Transcripts with CDS annotation: "
        f"{len(cds_by_transcript):,}"
    )

    print("CDS phase counts:")
    print(phase_counts)

    return cds_by_transcript, cds_starts


# ============================================================
# FUNCTION 2
# COUNT CODON POSITIONS
# ============================================================

def count_codon_positions(first_position, length):

    """
    Count codon positions 1, 2 and 3
    within a continuous coding segment.

    first_position:
        codon position of first base

    length:
        number of bases
    """

    if length <= 0:
        return np.array(
            [0, 0, 0],
            dtype=np.uint32
        )

    complete_codons, remainder = divmod(
        length,
        3
    )

    counts = np.array(
        [
            complete_codons,
            complete_codons,
            complete_codons
        ],
        dtype=np.uint32
    )

    for i in range(remainder):

        position = (
            first_position - 1 + i
        ) % 3

        counts[position] += np.uint32(1)

    return counts


# ============================================================
# FUNCTION 3
# ANNOTATE A GENOMIC RANGE
# ============================================================

def annotate_range(
    transcript_id,
    query_start,
    query_end,
    boundary_start,
    boundary_end,
    cds_by_transcript,
    cds_starts
):

    """
    Determine CDS overlap for a small genomic range.

    Returns:

        total_cds
        boundary_cds
        total_codon_counts
        boundary_codon_counts
    """

    total_cds = 0
    boundary_cds = 0

    total_codon = np.zeros(
        3,
        dtype=np.uint32
    )

    boundary_codon = np.zeros(
        3,
        dtype=np.uint32
    )

    intervals = cds_by_transcript.get(
        transcript_id
    )

    if intervals is None:

        return (
            0,
            0,
            np.zeros(3, dtype=np.uint32),
            np.zeros(3, dtype=np.uint32)
        )

    starts = cds_starts[transcript_id]

    # --------------------------------------------------------
    # Find first possible overlapping CDS interval
    # --------------------------------------------------------

    index = bisect.bisect_left(
        starts,
        query_start
    )

    if index > 0:
        index -= 1

    # --------------------------------------------------------
    # Examine overlapping CDS intervals
    # --------------------------------------------------------

    while index < len(intervals):

        cds_start, cds_end, strand, phase = (
            intervals[index]
        )

        # CDS begins after query region
        if cds_start > query_end:
            break

        # CDS ends before query region
        if cds_end < query_start:

            index += 1
            continue

        # ----------------------------------------------------
        # Determine overlap
        # ----------------------------------------------------

        overlap_start = max(
            query_start,
            cds_start
        )

        overlap_end = min(
            query_end,
            cds_end
        )

        if overlap_start <= overlap_end:

            length = (
                overlap_end -
                overlap_start +
                1
            )

            total_cds += length

            # ------------------------------------------------
            # Calculate codon positions
            # ------------------------------------------------

            if phase is not None:

                if strand == "+":

                    offset = (
                        overlap_start -
                        cds_start
                    )

                else:

                    offset = (
                        cds_end -
                        overlap_end
                    )

                first_position = (
                    (offset - phase) % 3
                ) + 1

                counts = count_codon_positions(
                    first_position,
                    length
                )

                total_codon += counts

            # ------------------------------------------------
            # Boundary overlap
            # ------------------------------------------------

            b_start = max(
                overlap_start,
                boundary_start
            )

            b_end = min(
                overlap_end,
                boundary_end
            )

            if b_start <= b_end:

                boundary_length = (
                    b_end -
                    b_start +
                    1
                )

                boundary_cds += (
                    boundary_length
                )

                if phase is not None:

                    if strand == "+":

                        boundary_offset = (
                            b_start -
                            cds_start
                        )

                    else:

                        boundary_offset = (
                            cds_end -
                            b_end
                        )

                    boundary_first_position = (
                        (
                            boundary_offset -
                            phase
                        ) % 3
                    ) + 1

                    boundary_counts = (
                        count_codon_positions(
                            boundary_first_position,
                            boundary_length
                        )
                    )

                    boundary_codon += (
                        boundary_counts
                    )

        index += 1

    return (
        total_cds,
        boundary_cds,
        total_codon,
        boundary_codon
    )


# ============================================================
# MAIN
# ============================================================

def main():

    # ========================================================
    # LOAD CDS
    # ========================================================

    cds_by_transcript, cds_starts = (
        load_cds_annotation()
    )

    # ========================================================
    # LOAD MASTER TABLE
    # ========================================================

    print()
    print("=" * 70)
    print("LOADING MASTER GENE-JUNCTION TABLE")
    print("=" * 70)

    master = pd.read_csv(
        MASTER,
        sep="\t",
        dtype=str
    )

    # Convert required numeric columns
    master["architecture_fraction"] = pd.to_numeric(
        master["architecture_fraction"],
        errors="coerce"
    )

    master["transcripts"] = pd.to_numeric(
        master["transcripts"],
        errors="coerce"
    ).fillna(0).astype(np.uint32)

    print(
        f"Master rows: {len(master):,}"
    )

    # ========================================================
    # MASTER LOOKUP
    # ========================================================

    print()
    print("Building master lookup...")

    master_lookup = {}

    for i, row in enumerate(
        master[
            ["gene_id", "junction_id"]
        ].itertuples(index=False)
    ):

        master_lookup[
            (
                row.gene_id,
                row.junction_id
            )
        ] = i

    n_junctions = len(master)

    print(
        f"Unique gene-junction records: "
        f"{n_junctions:,}"
    )

    # ========================================================
    # ARRAYS
    # ========================================================

    support_seen = np.zeros(
        n_junctions,
        dtype=np.uint32
    )

    # --------------------------------------------------------
    # CDS base counts
    # --------------------------------------------------------

    donor_flank_cds = np.zeros(
        n_junctions,
        dtype=np.uint32
    )

    donor_boundary_cds = np.zeros(
        n_junctions,
        dtype=np.uint32
    )

    acceptor_flank_cds = np.zeros(
        n_junctions,
        dtype=np.uint32
    )

    acceptor_boundary_cds = np.zeros(
        n_junctions,
        dtype=np.uint32
    )

    # --------------------------------------------------------
    # Codon counts
    # --------------------------------------------------------

    donor_flank_codon = np.zeros(
        (n_junctions, 3),
        dtype=np.uint32
    )

    donor_boundary_codon = np.zeros(
        (n_junctions, 3),
        dtype=np.uint32
    )

    acceptor_flank_codon = np.zeros(
        (n_junctions, 3),
        dtype=np.uint32
    )

    acceptor_boundary_codon = np.zeros(
        (n_junctions, 3),
        dtype=np.uint32
    )

    # ========================================================
    # PROCESS SPLICE JUNCTIONS
    # ========================================================

    print()
    print("=" * 70)
    print("PROCESSING TRANSCRIPT-SPECIFIC JUNCTIONS")
    print("=" * 70)

    rows_processed = 0
    missing_master = 0

    reader = pd.read_csv(
        SPLICE_JUNCTIONS,
        sep="\t",
        dtype={
            "gene_id": str,
            "transcript_id": str,
            "chromosome": str,
            "strand": str
        },
        chunksize=CHUNK_SIZE
    )

    for chunk_number, chunk in enumerate(
        reader,
        start=1
    ):

        for row in chunk.itertuples(
            index=False
        ):

            gene_id = row.gene_id
            transcript_id = row.transcript_id
            chromosome = row.chromosome
            strand = row.strand

            donor = int(row.donor)
            acceptor = int(row.acceptor)

            # ------------------------------------------------
            # Build junction ID
            # ------------------------------------------------

            left = min(
                donor,
                acceptor
            )

            right = max(
                donor,
                acceptor
            )

            junction_id = (
                f"{chromosome}:"
                f"{left}-{right}:"
                f"{strand}"
            )

            key = (
                gene_id,
                junction_id
            )

            idx = master_lookup.get(key)

            if idx is None:

                missing_master += 1
                continue

            support_seen[idx] += np.uint32(1)

            # =================================================
            # DONOR EXONIC REGION
            # =================================================

            if strand == "+":

                # -20 ... -1
                donor_start = donor - 20
                donor_end = donor - 1

                # -2 ... -1
                donor_boundary_start = donor - 2
                donor_boundary_end = donor - 1

            else:

                # -20 ... -1
                donor_start = donor + 1
                donor_end = donor + 20

                # -2 ... -1
                donor_boundary_start = donor + 1
                donor_boundary_end = donor + 2

            donor_query_start = min(
                donor_start,
                donor_end
            )

            donor_query_end = max(
                donor_start,
                donor_end
            )

            donor_boundary_query_start = min(
                donor_boundary_start,
                donor_boundary_end
            )

            donor_boundary_query_end = max(
                donor_boundary_start,
                donor_boundary_end
            )

            (
                donor_total_cds,
                donor_boundary_cds_i,
                donor_total_codon,
                donor_boundary_codon_i
            ) = annotate_range(
                transcript_id,
                donor_query_start,
                donor_query_end,
                donor_boundary_query_start,
                donor_boundary_query_end,
                cds_by_transcript,
                cds_starts
            )

            donor_flank_cds_i = (
                donor_total_cds -
                donor_boundary_cds_i
            )

            donor_flank_codon_i = (
                donor_total_codon -
                donor_boundary_codon_i
            )

            donor_flank_cds[idx] += np.uint32(
                donor_flank_cds_i
            )

            donor_boundary_cds[idx] += np.uint32(
                donor_boundary_cds_i
            )

            donor_flank_codon[idx] += (
                donor_flank_codon_i.astype(
                    np.uint32
                )
            )

            donor_boundary_codon[idx] += (
                donor_boundary_codon_i.astype(
                    np.uint32
                )
            )

            # =================================================
            # ACCEPTOR EXONIC REGION
            # =================================================

            if strand == "+":

                # +1 ... +20
                acceptor_start = acceptor + 1
                acceptor_end = acceptor + 20

                # +1 ... +2
                acceptor_boundary_start = acceptor + 1
                acceptor_boundary_end = acceptor + 2

            else:

                # +1 ... +20
                acceptor_start = acceptor - 20
                acceptor_end = acceptor - 1

                # +1 ... +2
                acceptor_boundary_start = acceptor - 2
                acceptor_boundary_end = acceptor - 1

            acceptor_query_start = min(
                acceptor_start,
                acceptor_end
            )

            acceptor_query_end = max(
                acceptor_start,
                acceptor_end
            )

            acceptor_boundary_query_start = min(
                acceptor_boundary_start,
                acceptor_boundary_end
            )

            acceptor_boundary_query_end = max(
                acceptor_boundary_start,
                acceptor_boundary_end
            )

            (
                acceptor_total_cds,
                acceptor_boundary_cds_i,
                acceptor_total_codon,
                acceptor_boundary_codon_i
            ) = annotate_range(
                transcript_id,
                acceptor_query_start,
                acceptor_query_end,
                acceptor_boundary_query_start,
                acceptor_boundary_query_end,
                cds_by_transcript,
                cds_starts
            )

            acceptor_flank_cds_i = (
                acceptor_total_cds -
                acceptor_boundary_cds_i
            )

            acceptor_flank_codon_i = (
                acceptor_total_codon -
                acceptor_boundary_codon_i
            )

            acceptor_flank_cds[idx] += np.uint32(
                acceptor_flank_cds_i
            )

            acceptor_boundary_cds[idx] += np.uint32(
                acceptor_boundary_cds_i
            )

            acceptor_flank_codon[idx] += (
                acceptor_flank_codon_i.astype(
                    np.uint32
                )
            )

            acceptor_boundary_codon[idx] += (
                acceptor_boundary_codon_i.astype(
                    np.uint32
                )
            )

            rows_processed += 1

        print(
            f"Chunk {chunk_number}: "
            f"{rows_processed:,} "
            f"transcript-junctions processed"
        )

    # ========================================================
    # VALIDATION
    # ========================================================

    print()
    print("=" * 70)
    print("VALIDATION")
    print("=" * 70)

    expected_support = (
        master["transcripts"]
        .to_numpy(dtype=np.uint32)
    )

    support_mismatch = int(
        np.sum(
            support_seen != expected_support
        )
    )

    print(
        f"Transcript-junction rows processed: "
        f"{rows_processed:,}"
    )

    print(
        f"Missing master records: "
        f"{missing_master:,}"
    )

    print(
        f"Support-count mismatches: "
        f"{support_mismatch:,}"
    )

    if support_mismatch == 0:

        print(
            "SUCCESS: support counts match "
            "the master table."
        )

    else:

        print(
            "WARNING: support counts do not match."
        )

    # ========================================================
    # CREATE OUTPUT
    # ========================================================

    print()
    print("=" * 70)
    print("CREATING CODING-CONTEXT TABLE")
    print("=" * 70)

    result = master[
        [
            "gene_id",
            "gene_name",
            "junction_id",
            "architecture_fraction",
            "transcripts"
        ]
    ].copy()

    support = expected_support.astype(
        np.float64
    )

    # --------------------------------------------------------
    # CDS FRACTIONS
    # --------------------------------------------------------

    result[
        "donor_exonic_flank_cds_fraction"
    ] = np.divide(
        donor_flank_cds,
        support * 18,
        out=np.zeros(
            n_junctions,
            dtype=float
        ),
        where=support > 0
    )

    result[
        "donor_boundary_cds_fraction"
    ] = np.divide(
        donor_boundary_cds,
        support * 2,
        out=np.zeros(
            n_junctions,
            dtype=float
        ),
        where=support > 0
    )

    result[
        "acceptor_exonic_flank_cds_fraction"
    ] = np.divide(
        acceptor_flank_cds,
        support * 18,
        out=np.zeros(
            n_junctions,
            dtype=float
        ),
        where=support > 0
    )

    result[
        "acceptor_boundary_cds_fraction"
    ] = np.divide(
        acceptor_boundary_cds,
        support * 2,
        out=np.zeros(
            n_junctions,
            dtype=float
        ),
        where=support > 0
    )

    # ========================================================
    # CODON FRACTION FUNCTION
    # ========================================================

    def add_codon_fractions(
        dataframe,
        prefix,
        codon_counts,
        coding_counts
    ):

        denominator = (
            coding_counts.astype(
                np.float64
            )
        )

        for j in range(3):

            codon_number = j + 1

            dataframe[
                f"{prefix}_codon{codon_number}_fraction"
            ] = np.divide(
                codon_counts[:, j],
                denominator,
                out=np.zeros(
                    n_junctions,
                    dtype=float
                ),
                where=denominator > 0
            )

    # --------------------------------------------------------
    # Add codon fractions
    # --------------------------------------------------------

    add_codon_fractions(
        result,
        "donor_exonic_flank",
        donor_flank_codon,
        donor_flank_cds
    )

    add_codon_fractions(
        result,
        "donor_boundary",
        donor_boundary_codon,
        donor_boundary_cds
    )

    add_codon_fractions(
        result,
        "acceptor_exonic_flank",
        acceptor_flank_codon,
        acceptor_flank_cds
    )

    add_codon_fractions(
        result,
        "acceptor_boundary",
        acceptor_boundary_codon,
        acceptor_boundary_cds
    )

    # ========================================================
    # SAVE
    # ========================================================

    result.to_csv(
        OUTPUT,
        sep="\t",
        index=False
    )

    print()
    print("=" * 70)
    print("FINISHED")
    print("=" * 70)

    print(
        f"Output rows: {len(result):,}"
    )

    print(
        f"Output columns: {len(result.columns):,}"
    )

    print(
        f"Output file: {OUTPUT}"
    )

    print()
    print("First five rows:")

    print(
        result.head().to_string(
            index=False
        )
    )


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":

    main()
