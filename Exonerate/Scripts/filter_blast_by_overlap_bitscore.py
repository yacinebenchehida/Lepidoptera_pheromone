#!/usr/bin/env python3

import sys
from collections import defaultdict

def parse_blast_line(line):
    parts = line.strip().split("\t")
    chrom_info = parts[0]
    bit_score = float(parts[11])
    chrom, start, end, *_ = chrom_info.split("_")
    start, end = int(start), int(end)
    return chrom, start, end, bit_score, line.strip()

def intervals_overlap(a_start, a_end, b_start, b_end):
    return not (a_end < b_start or b_end < a_start)

def process_blast_file(blast_file):
    chrom_hits = defaultdict(list)

    with open(blast_file) as f:
        for line in f:
            chrom, start, end, bit_score, raw_line = parse_blast_line(line)
            chrom_hits[chrom].append((start, end, bit_score, raw_line))

    retained = []
    for chrom, entries in chrom_hits.items():
        entries.sort(key=lambda x: x[0])  # sort by start
        groups = []

        for entry in entries:
            placed = False
            for group in groups:
                if any(intervals_overlap(entry[0], entry[1], g[0], g[1]) for g in group):
                    group.append(entry)
                    placed = True
                    break
            if not placed:
                groups.append([entry])

        for group in groups:
            best = max(group, key=lambda x: x[2])  # highest bit score
            retained.append(best[3])

    return retained

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python filter_blast_by_overlap_bitscore.py <blast_output.tsv>")
        sys.exit(1)

    filtered = process_blast_file(sys.argv[1])
    for line in filtered:
        print(line)

