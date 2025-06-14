#!/usr/bin/env python

import sys  # Importing the sys module to access command-line arguments
import os  # Importing the os module to perform operating system related tasks
from Bio import SeqIO  # Importing SeqIO module from Biopython for sequence parsing

#################################################################
# Function to write chunks of sequences to separate FASTA files #
#################################################################
def write_chunks(input_file):
    # Get the directory of the input file
    input_dir = os.path.dirname(input_file)

    # Open the input FASTA file for reading
    with open(input_file, "r") as f:
        # Parse the input FASTA file and store the records in a list
        records = list(SeqIO.parse(f, "fasta"))

        # Calculate the number of chunks based on the number of records
        num_chunks = len(records) // 10 + (1 if len(records) % 10 != 0 else 0)

        # Iterate over each chunk
        for i in range(num_chunks):
            # Calculate the start and end index of the current chunk
            chunk_start = i * 10
            chunk_end = min((i + 1) * 10, len(records))
            # Get the records for the current chunk
            chunk_records = records[chunk_start:chunk_end]

            # Construct the output file path for the current chunk
            output_file = os.path.join(input_dir, f"output_chunk_{i + 1}.fasta")

            # Write the chunk of sequences to a new FASTA file
            with open(output_file, "w") as out_f:
                SeqIO.write(chunk_records, out_f, "fasta")

# Check if the script is being run as the main program
if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 2:
        # Print usage message and exit with status 1 if incorrect number of arguments
        print("Usage: python script.py input.fasta")
        sys.exit(1)
    # Get input file path from command-line arguments
    input_file = sys.argv[1]
    # Call the write_chunks function with provided input file
    write_chunks(input_file)
