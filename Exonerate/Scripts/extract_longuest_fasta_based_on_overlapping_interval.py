#!/usr/bin/env python3
import sys
import re
from collections import defaultdict
from Bio import SeqIO

def parse_header(header):
    # Extract: chrom, start, end from the header
    match = re.match(r"^(\S+?)_(\d+)_(\d+)", header)
    if match:
        chrom = match.group(1)
        start = int(match.group(2))
        end = int(match.group(3))
        return chrom, start, end
    else:
        raise ValueError(f"Header format not recognized: {header}")

def intervals_overlap(start1, end1, start2, end2):
    return max(start1, start2) <= min(end1, end2)

def main(fasta_file):
    chrom_intervals = defaultdict(list)

    for record in SeqIO.parse(fasta_file, "fasta"):
        chrom, start, end = parse_header(record.id)
        chrom_intervals[chrom].append((start, end, record))

    for chrom in chrom_intervals:
        intervals = chrom_intervals[chrom]
        intervals.sort(key=lambda x: x[0])  # sort by start position
        selected = []

        while intervals:
            base_start, base_end, base_record = intervals.pop(0)
            overlap_group = [(base_start, base_end, base_record)]

            non_overlapping = []
            for i_start, i_end, i_record in intervals:
                if intervals_overlap(base_start, base_end, i_start, i_end):
                    overlap_group.append((i_start, i_end, i_record))
                else:
                    non_overlapping.append((i_start, i_end, i_record))

            # choose longest sequence in overlap group
            longest = max(overlap_group, key=lambda x: len(x[2].seq))
            selected.append(longest[2])
            intervals = non_overlapping

        for rec in selected:
            print(f">{rec.id} {rec.description[len(rec.id):].strip()}")
            print(str(rec.seq))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.fasta>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])

