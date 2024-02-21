#!/usr/bin/env python

import sys
import os
from Bio import SeqIO

def write_chunks(input_file):
	input_dir = os.path.dirname(input_file)

	# Open the input FASTA file for reading
	with open(input_file, "r") as f:
		records = list(SeqIO.parse(f, "fasta"))

		num_chunks = len(records) // 100 + (1 if len(records) % 100 != 0 else 0)

		for i in range(num_chunks):
			chunk_start = i * 100
			chunk_end = min((i + 1) * 100, len(records))
			chunk_records = records[chunk_start:chunk_end]

			output_file = os.path.join(input_dir, f"output_chunk_{i + 1}.fasta")

			# Write the chunk of sequences to a new FASTA file
			with open(output_file, "w") as out_f:
				SeqIO.write(chunk_records, out_f, "fasta")

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage: python script.py input.fasta")
		sys.exit(1)
	input_file = sys.argv[1]
	write_chunks(input_file)
