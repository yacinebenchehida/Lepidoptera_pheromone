#!/usr/bin/env python3
import sys
import random
from Bio import SeqIO

def main(fasta_file):
    seq_to_records = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_str = str(record.seq)
        if seq_str not in seq_to_records:
            seq_to_records[seq_str] = []
        seq_to_records[seq_str].append(record)

    for seq_str, records in seq_to_records.items():
        chosen_record = random.choice(records)
        print(f">{chosen_record.id} {chosen_record.description[len(chosen_record.id):].strip()}")
        print(seq_str)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <input.fasta>", file=sys.stderr)
        sys.exit(1)
    main(sys.argv[1])









