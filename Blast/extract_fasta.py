#!/usr/bin/env python3
import sys
from Bio import SeqIO

#############################################################################################################################
# Function that fishes specific sequences when their headers contain a specified pattern (taken from a larger fasta file)    #
#############################################################################################################################
def extract_sequences(id_file, fasta_file, output_file):
    """
    Extracts sequences from a FASTA file if the header contains any pattern from the ID file.

    Args:
        id_file (str): Path to the file containing the ID patterns.
        fasta_file (str): Path to the input FASTA file.
        output_file (str): Path to the output file to write the extracted sequences.

    Returns:
        None
    """
    # Read IDs (patterns) from the file
    with open(id_file, 'r') as f:
        ids = {line.strip() for line in f if line.strip()}  # Remove empty lines and strip whitespace

    # Extract sequences whose headers contain any of the patterns
    with open(output_file, 'w') as out_f:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            # Check if any ID pattern is a substring of the header
            if any(id_pattern in record.description for id_pattern in ids):
                SeqIO.write(record, out_f, 'fasta')

if __name__ == "__main__":
    # Check if the correct number of command-line arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python extract_fasta.py <id_file> <fasta_file> <output_file>")
        sys.exit(1)
    
    # Parse command-line arguments
    id_file = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]

    # Call the function to extract sequences
    extract_sequences(id_file, fasta_file, output_file)
