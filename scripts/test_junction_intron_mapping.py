import pandas as pd


# ==========================================================
# 1. READ THE TWO TABLES
# ==========================================================

events = pd.read_csv(
    "data/processed/validated_exon_skipping_events.tsv",
    sep="\t"
)

sites = pd.read_csv(
    "data/processed/unique_splice_sites.tsv",
    sep="\t"
)


# ==========================================================
# 2. CREATE A FAST LOOKUP OF UNIQUE SPLICE SITES
# ==========================================================

site_ids = set(
    sites["intron_id"]
)


# ==========================================================
# 3. TEST FIRST 20 UNIQUE SKIPPING JUNCTIONS
# ==========================================================

rows = []


unique_events = (
    events
    .drop_duplicates("skipping_junction")
    .head(20)
)


for _, row in unique_events.iterrows():

    junction = row["skipping_junction"]

    chromosome = row["chromosome"]

    strand = row["strand"]


    # ------------------------------------------------------
    # Extract donor and acceptor from:
    #
    # chrX:100633405-100630866:-
    # ------------------------------------------------------

    coordinate_part = (
        junction
        .split(":")[1]
    )

    donor = int(
        coordinate_part
        .split("-")[0]
    )

    acceptor = int(
        coordinate_part
        .split("-")[1]
    )


    # ------------------------------------------------------
    # Convert junction coordinates
    # into intron coordinates
    # ------------------------------------------------------

    if strand == "-":

        intron_start = acceptor + 1
        intron_end = donor - 1

    else:

        intron_start = donor + 1
        intron_end = acceptor - 1


    # ------------------------------------------------------
    # Build intron ID
    # ------------------------------------------------------

    intron_id = (
        f"{chromosome}:"
        f"{intron_start}-"
        f"{intron_end}:"
        f"{strand}"
    )


    # ------------------------------------------------------
    # Check whether it exists
    # in unique_splice_sites.tsv
    # ------------------------------------------------------

    found = (
        intron_id in site_ids
    )


    rows.append([

        junction,

        intron_id,

        found

    ])


# ==========================================================
# 4. DISPLAY RESULTS
# ==========================================================

result = pd.DataFrame(

    rows,

    columns=[

        "skipping_junction",

        "mapped_intron_id",

        "found_in_unique_sites"

    ]

)


print()

print(
    result.to_string(index=False)
)


print()

print(
    "Tested:",
    len(result)
)

print(
    "Found:",
    result["found_in_unique_sites"].sum()
)

print(
    "Not found:",
    (~result["found_in_unique_sites"]).sum()
)
