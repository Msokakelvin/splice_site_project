import pandas as pd


event_file = "data/processed/strict_exon_skipping_events.tsv"
junction_file = "data/processed/splice_junctions.tsv"
output_file = "data/processed/validated_exon_skipping_events.tsv"


print("Reading junction table...")


# ==========================================================
# 1. READ JUNCTION TABLE
# ==========================================================

junctions = pd.read_csv(
    junction_file,
    sep="\t",
    usecols=[
        "gene_id",
        "transcript_id",
        "chromosome",
        "strand",
        "donor",
        "acceptor"
    ]
)

print("Junction rows:", len(junctions))


# ==========================================================
# 2. TRANSCRIPT INFORMATION
# ==========================================================

transcript_info = (
    junctions[
        [
            "transcript_id",
            "chromosome",
            "strand"
        ]
    ]
    .drop_duplicates()
)


transcript_strand = (
    transcript_info
    .set_index("transcript_id")["strand"]
    .to_dict()
)


transcript_chromosome = (
    transcript_info
    .set_index("transcript_id")["chromosome"]
    .to_dict()
)


print(
    "Transcript strand entries:",
    len(transcript_strand)
)

print(
    "Transcript chromosome entries:",
    len(transcript_chromosome)
)


# ==========================================================
# 3. CREATE FAST JUNCTION LOOKUP
# ==========================================================

junction_lookup = set(
    zip(
        junctions["transcript_id"],
        junctions["donor"],
        junctions["acceptor"]
    )
)


print(
    "Unique transcript junctions:",
    len(junction_lookup)
)


# ==========================================================
# 4. READ EXON-SKIPPING EVENTS
# ==========================================================

print("Reading exon-skipping events...")


events = pd.read_csv(
    event_file,
    sep="\t"
)

events = events.drop_duplicates()


print(
    "Events to validate:",
    len(events)
)


# ==========================================================
# 5. VALIDATE EVENTS
# ==========================================================

results = []


for _, row in events.iterrows():

    transcript_included = row[
        "transcript_with_exon"
    ]

    transcript_skipping = row[
        "transcript_skipping_exon"
    ]


    strand = transcript_strand.get(
        transcript_included
    )

    chromosome = transcript_chromosome.get(
        transcript_included
    )


    if strand is None:
        continue

    if chromosome is None:
        continue


    upstream_start = int(
        row["upstream_exon_start"]
    )

    upstream_end = int(
        row["upstream_exon_end"]
    )

    skipped_start = int(
        row["skipped_exon_start"]
    )

    skipped_end = int(
        row["skipped_exon_end"]
    )

    downstream_start = int(
        row["downstream_exon_start"]
    )

    downstream_end = int(
        row["downstream_exon_end"]
    )


    # ======================================================
    # POSITIVE STRAND
    # ======================================================

    if strand == "+":

        junction_1 = (
            transcript_included,
            upstream_end,
            skipped_start
        )

        junction_2 = (
            transcript_included,
            skipped_end,
            downstream_start
        )

        junction_3 = (
            transcript_skipping,
            upstream_end,
            downstream_start
        )


    # ======================================================
    # NEGATIVE STRAND
    # ======================================================

    else:

        junction_1 = (
            transcript_included,
            upstream_start,
            skipped_end
        )

        junction_2 = (
            transcript_included,
            skipped_start,
            downstream_end
        )

        junction_3 = (
            transcript_skipping,
            upstream_start,
            downstream_end
        )


    # ======================================================
    # CHECK JUNCTIONS
    # ======================================================

    has_junction_1 = (
        junction_1 in junction_lookup
    )

    has_junction_2 = (
        junction_2 in junction_lookup
    )

    has_junction_3 = (
        junction_3 in junction_lookup
    )


    fully_validated = (
        has_junction_1
        and has_junction_2
        and has_junction_3
    )


    # ======================================================
    # CREATE ACTUAL SKIPPING JUNCTION ID
    # ======================================================

    skipping_junction_id = (
        str(chromosome)
        + ":"
        + str(int(junction_3[1]))
        + "-"
        + str(int(junction_3[2]))
        + ":"
        + str(strand)
    )


    # ======================================================
    # SAVE RESULT
    # ======================================================

    results.append([

        row["gene_id"],

        transcript_included,

        transcript_skipping,

        chromosome,

        strand,

        skipping_junction_id,

        has_junction_1,

        has_junction_2,

        has_junction_3,

        fully_validated

    ])


# ==========================================================
# 6. CREATE OUTPUT TABLE
# ==========================================================

result = pd.DataFrame(
    results,
    columns=[

        "gene_id",

        "transcript_with_exon",

        "transcript_skipping_exon",

        "chromosome",

        "strand",

        "skipping_junction",

        "included_junction_1",

        "included_junction_2",

        "skipping_junction_validated",

        "fully_validated"

    ]
)


# ==========================================================
# 7. SAVE OUTPUT
# ==========================================================

result.to_csv(
    output_file,
    sep="\t",
    index=False
)


# ==========================================================
# 8. SUMMARY
# ==========================================================

print()
print("====================================")
print("VALIDATION COMPLETE")
print("====================================")

print(
    "Events checked:",
    len(result)
)

print(
    "Fully validated:",
    result["fully_validated"].sum()
)

print(
    "Not fully validated:",
    (~result["fully_validated"]).sum()
)

print(
    "Validation rate:",
    round(
        result["fully_validated"].mean() * 100,
        2
    ),
    "%"
)

print(
    "Unique skipping junctions:",
    result["skipping_junction"].nunique()
)

print(
    "Unique gene + skipping junction:",
    result[
        ["gene_id", "skipping_junction"]
    ]
    .drop_duplicates()
    .shape[0]
)

print(
    "Unique genes:",
    result["gene_id"].nunique()
)

print(
    "Output:",
    output_file
)
