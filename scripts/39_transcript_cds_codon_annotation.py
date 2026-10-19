#!/usr/bin/env python3
"""
Script 39 – Transcript-aware CDS / codon-phase annotation for splice junctions.
"""

import gzip
import os
import sys
import bisect
from collections import defaultdict
import pandas as pd

# ── paths ──────────────────────────────────────────────────────────────
JUNCTION_MAP = "data/processed/junction_transcript_map.tsv"
OUTPUT_FILE  = "data/processed/junction_position_cds_summary.tsv"

# >>>  CORRECTED PATH  <<<
GENCODE_GTF  = "data/reference/gencode.v50.annotation.gtf.gz"

# ── positions to annotate (exon side only) ────────────────────────────
DONOR_EXON_POS    = list(range(-20, 0))   # −20 … −1  (upstream exon)
ACCEPTOR_EXON_POS = list(range(1, 21))    #  +1 … +20  (downstream exon)

# ── helpers ────────────────────────────────────────────────────────────
def load_cds_from_gtf(gtf_path):
    """Return dict  transcript_id → [(chrom, start, end, strand), …]"""
    print(f"Loading CDS intervals from {gtf_path}...")
    cds = defaultdict(list)

    if not os.path.exists(gtf_path):
        sys.exit(f"ERROR: GTF not found at {gtf_path}")

    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 9 or cols[2] != "CDS":
                continue

            chrom, start, end = cols[0], int(cols[3]), int(cols[4])
            strand = cols[6]

            tid = None
            for attr in cols[8].split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id"):
                    tid = attr.split('"')[1]
                    break

            if tid:
                cds[tid].append((chrom, start, end, strand))

    print(f"  Found CDS for {len(cds):,} transcripts.")
    return cds

def build_lookup(intervals):
    """Build bisect-ready arrays from a transcript's CDS intervals."""
    chrom = intervals[0][0]
    strand = intervals[0][3]

    intervals.sort(key=lambda x: x[1])

    starts = [iv[1] for iv in intervals]
    ends   = [iv[2] for iv in intervals]

    reading_intervals = intervals if strand == "+" else intervals[::-1]

    cum_in_reading_order = []
    current_cum = 0
    for iv in reading_intervals:
        cum_in_reading_order.append(current_cum)
        current_cum += (iv[2] - iv[1] + 1)

    if strand == "-":
        cum_in_reading_order = cum_in_reading_order[::-1]

    return chrom, strand, starts, ends, cum_in_reading_order

def query(pos, chrom, strand, starts, ends, cum_offsets):
    """Return (in_cds: bool, codon_pos: int|None)."""
    idx = bisect.bisect_right(starts, pos) - 1
    if idx >= 0 and pos <= ends[idx]:
        offset = (pos - starts[idx]) if strand == "+" else (ends[idx] - pos)
        cds_pos = cum_offsets[idx] + offset
        codon  = (cds_pos % 3) + 1
        return True, codon
    return False, None

# ── main ───────────────────────────────────────────────────────────────
def main():
    print("=== Script 39: Transcript-aware CDS / codon-phase annotation ===\n")

    cds_by_tid = load_cds_from_gtf(GENCODE_GTF)

    print(f"\nLoading {JUNCTION_MAP}...")
    jmap = pd.read_csv(JUNCTION_MAP, sep="\t", dtype=str)
    print(f"  Loaded {len(jmap):,} junction–transcript pairs.")

    print("Parsing junction coordinates...")
    jmap["j_chrom"] = jmap["junction_id"].str.split(":").str[0]
    coords = jmap["junction_id"].str.split(":").str[1].str.split("-")
    jmap["j_start"] = coords.str[0].astype(int)
    jmap["j_end"] = coords.str[1].astype(int)
    jmap["j_strand"] = jmap["junction_id"].str.split(":").str[2]

    counts = defaultdict(lambda: [0, 0, 0, 0, 0])

    grouped = jmap.groupby("transcript_id")
    n_groups = len(grouped)
    done_transcripts = 0
    done_junctions = 0

    print(f"\nAnnotating (grouped by {n_groups:,} transcripts)...")

    for tid, grp in grouped:
        intervals = cds_by_tid.get(tid)
        if intervals is None:
            tid_nover = tid.split(".")[0]
            intervals = cds_by_tid.get(tid_nover)

        has_cds = False
        if intervals:
            chrom, strand, starts, ends, cum = build_lookup(intervals)
            has_cds = True

        for jid, j_chrom, j_start, j_end, j_strand in grp[["junction_id", "j_chrom", "j_start", "j_end", "j_strand"]].itertuples(index=False, name=None):
            if j_strand == "+":
                donor_g, acc_g = j_start, j_end
            else:
                donor_g, acc_g = j_end, j_start

            for p in DONOR_EXON_POS:
                gpos = donor_g + p if j_strand == "+" else donor_g - p
                row = counts[(jid, "donor", p)]
                row[0] += 1
                if has_cds and j_chrom == chrom:
                    in_cds, cp = query(gpos, j_chrom, strand, starts, ends, cum)
                    if in_cds:
                        row[1] += 1
                        row[1 + cp] += 1

            for p in ACCEPTOR_EXON_POS:
                gpos = acc_g + p if j_strand == "+" else acc_g - p
                row = counts[(jid, "acceptor", p)]
                row[0] += 1
                if has_cds and j_chrom == chrom:
                    in_cds, cp = query(gpos, j_chrom, strand, starts, ends, cum)
                    if in_cds:
                        row[1] += 1
                        row[1 + cp] += 1

            done_junctions += 1

        done_transcripts += 1
        if done_transcripts % 10000 == 0:
            print(f"  {done_transcripts:>8,} transcripts processed ({done_transcripts / n_groups * 100:.1f}%)")

    print(f"\n  Finished: {done_junctions:,} junctions annotated across {done_transcripts:,} transcripts.")

    print(f"\nWriting {OUTPUT_FILE}...")
    with open(OUTPUT_FILE, "w") as fh:
        fh.write("junction_id\tsite\tposition\ttotal\tcds\tnon_cds\tcds_frac\tcodon1\tcodon2\tcodon3\n")
        for (jid, site, pos), v in sorted(counts.items()):
            tot, cds, c1, c2, c3 = v
            nc = tot - cds
            frac = cds / tot if tot else 0.0
            fh.write(f"{jid}\t{site}\t{pos}\t{tot}\t{cds}\t{nc}\t{frac:.4f}\t{c1}\t{c2}\t{c3}\n")

    print(f"  {len(counts):,} rows written.")
    print("\n✅ Script 39 complete.")

if __name__ == "__main__":
    main()
