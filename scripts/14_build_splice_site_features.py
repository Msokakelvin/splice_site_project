import pandas as pd


# ==========================================================
# FILES
# ==========================================================

site_file = (
    "data/processed/unique_splice_sites.tsv"
)

junction_file = (
    "data/processed/splice_junctions.tsv"
)

skipping_file = (
    "data/processed/validated_exon_skipping_events.tsv"
)

alternative_file = (
    "data/processed/alternative_donor_acceptor_events.tsv"
)

output_file = (
    "data/processed/splice_site_features.tsv"
)


# ==========================================================
# 1. READ UNIQUE SPLICE SITES
# ==========================================================

print("Reading unique splice sites...")

sites = pd.read_csv(
    site_file,
    sep="\t"
)

print(
    "Unique splice sites:",
    len(sites)
)


# ==========================================================
# 2. READ SPLICE JUNCTIONS
# ==========================================================

print()
print("Reading splice junctions...")

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

print(
    "Junction rows:",
    len(junctions)
)


# ==========================================================
# 3. CREATE JUNCTION ID
# ==========================================================

junctions["junction_id"] = (
    junctions["chromosome"].astype(str)
    + ":"
    + junctions["donor"].astype(int).astype(str)
    + "-"
    + junctions["acceptor"].astype(int).astype(str)
    + ":"
    + junctions["strand"].astype(str)
)


# ==========================================================
# 4. EXON-SKIPPING JUNCTIONS
# ==========================================================

print()
print("Reading exon-skipping events...")

skipping = pd.read_csv(
    skipping_file,
    sep="\t"
)

print(
    "Exon-skipping rows:",
    len(skipping)
)


skipping_junctions = set(
    skipping[
        "skipping_junction"
    ]
    .dropna()
    .unique()
)


print(
    "Unique exon-skipping junctions:",
    len(skipping_junctions)
)


# ==========================================================
# 5. ALTERNATIVE DONOR GROUPS
# ==========================================================

print()
print(
    "Reading alternative donor/acceptor groups..."
)

alternative = pd.read_csv(
    alternative_file,
    sep="\t"
)


# ----------------------------------------------------------
# Alternative donor
#
# Same acceptor + multiple donors
# ----------------------------------------------------------

alternative_donor = alternative[
    alternative["event_type"]
    == "alternative_5prime_donor"
].copy()


alternative_donor_acceptors = set(
    alternative_donor[
        "acceptor"
    ]
    .dropna()
    .astype(int)
)


print(
    "Alternative donor groups:",
    len(alternative_donor_acceptors)
)


# ----------------------------------------------------------
# Alternative acceptor
#
# Same donor + multiple acceptors
# ----------------------------------------------------------

alternative_acceptor = alternative[
    alternative["event_type"]
    == "alternative_3prime_acceptor"
].copy()


alternative_acceptor_donors = set(
    alternative_acceptor[
        "donor"
    ]
    .dropna()
    .astype(int)
)


print(
    "Alternative acceptor groups:",
    len(alternative_acceptor_donors)
)


# ==========================================================
# 6. IDENTIFY ACTUAL JUNCTIONS INVOLVED
# ==========================================================

# ----------------------------------------------------------
# Alternative donor:
#
# A junction participates if its acceptor coordinate
# belongs to an alternative-donor group.
# ----------------------------------------------------------

alternative_donor_junctions = set(
    junctions.loc[
        junctions["acceptor"].isin(
            alternative_donor_acceptors
        ),
        "junction_id"
    ]
)


# ----------------------------------------------------------
# Alternative acceptor:
#
# A junction participates if its donor coordinate
# belongs to an alternative-acceptor group.
# ----------------------------------------------------------

alternative_acceptor_junctions = set(
    junctions.loc[
        junctions["donor"].isin(
            alternative_acceptor_donors
        ),
        "junction_id"
    ]
)


print()
print(
    "Junctions associated with alternative donor groups:",
    len(alternative_donor_junctions)
)


print(
    "Junctions associated with alternative acceptor groups:",
    len(alternative_acceptor_junctions)
)


# ==========================================================
# 7. CONVERT SPLICE-SITE IDS
# ==========================================================

# unique_splice_sites.tsv uses intron IDs:
#
# chromosome:intron_start-intron_end:strand
#
# while junctions use:
#
# chromosome:donor-acceptor:strand
#
# Convert junction coordinates to intron coordinates.


def junction_to_intron(row):

    donor = int(
        row["donor"]
    )

    acceptor = int(
        row["acceptor"]
    )

    strand = row["strand"]


    if strand == "-":

        intron_start = acceptor + 1

        intron_end = donor - 1

    else:

        intron_start = donor + 1

        intron_end = acceptor - 1


    return (
        f"{row['chromosome']}:"
        f"{intron_start}-"
        f"{intron_end}:"
        f"{strand}"
    )


junctions["intron_id"] = (
    junctions.apply(
        junction_to_intron,
        axis=1
    )
)


# ==========================================================
# 8. CREATE JUNCTION → INTRON LOOKUP
# ==========================================================

junction_to_intron_lookup = dict(
    zip(
        junctions["junction_id"],
        junctions["intron_id"]
    )
)


# ==========================================================
# 9. CONVERT EVENT JUNCTIONS TO INTRON IDS
# ==========================================================

skipping_introns = set(
    junction_to_intron_lookup[j]
    for j in skipping_junctions
    if j in junction_to_intron_lookup
)


alternative_donor_introns = set(
    junction_to_intron_lookup[j]
    for j in alternative_donor_junctions
    if j in junction_to_intron_lookup
)


alternative_acceptor_introns = set(
    junction_to_intron_lookup[j]
    for j in alternative_acceptor_junctions
    if j in junction_to_intron_lookup
)


print()
print(
    "Exon-skipping introns:",
    len(skipping_introns)
)


print(
    "Alternative-donor introns:",
    len(alternative_donor_introns)
)


print(
    "Alternative-acceptor introns:",
    len(alternative_acceptor_introns)
)


# ==========================================================
# 10. ADD STRUCTURAL FLAGS TO UNIQUE SPLICE SITES
# ==========================================================

sites["exon_skipping"] = (
    sites["intron_id"]
    .isin(skipping_introns)
)


sites["alternative_donor"] = (
    sites["intron_id"]
    .isin(alternative_donor_introns)
)


sites["alternative_acceptor"] = (
    sites["intron_id"]
    .isin(alternative_acceptor_introns)
)


# ==========================================================
# 11. CREATE STRUCTURAL EVIDENCE CLASS
# ==========================================================

sites["structural_evidence"] = "none"


sites.loc[
    sites["exon_skipping"],
    "structural_evidence"
] = "exon_skipping"


sites.loc[
    sites["alternative_donor"],
    "structural_evidence"
] = "alternative_donor"


sites.loc[
    sites["alternative_acceptor"],
    "structural_evidence"
] = "alternative_acceptor"


# Multiple evidence types

multiple = (
    sites[
        [
            "exon_skipping",
            "alternative_donor",
            "alternative_acceptor"
        ]
    ]
    .sum(axis=1)
    > 1
)


sites.loc[
    multiple,
    "structural_evidence"
] = "multiple_event_types"


# ==========================================================
# 12. SAVE FINAL FEATURE TABLE
# ==========================================================

sites.to_csv(
    output_file,
    sep="\t",
    index=False
)


# ==========================================================
# 13. SUMMARY
# ==========================================================

print()
print(
    "===================================="
)

print(
    "FEATURE TABLE COMPLETE"
)

print(
    "===================================="
)


print(
    "Total splice sites:",
    len(sites)
)


print(
    "Exon-skipping sites:",
    sites["exon_skipping"].sum()
)


print(
    "Alternative donor sites:",
    sites["alternative_donor"].sum()
)


print(
    "Alternative acceptor sites:",
    sites["alternative_acceptor"].sum()
)


print()
print(
    "Structural evidence:"
)


print(
    sites[
        "structural_evidence"
    ]
    .value_counts()
)


print()
print(
    "Output:",
    output_file
)
