from Bio import SeqIO
import sys

def filter_fasta_by_length(fasta_file, length_threshold):
    for record in SeqIO.parse(fasta_file, "fasta"):
        if len(record.seq) >= length_threshold:
            print(f">{record.id}")
            print(record.seq)

# Get inputs from the command line
if len(sys.argv) != 3:
    print("Usage: python3 ./filter_fasta_by_length.py <fasta_file> <length_threshold>")
    sys.exit(1)

fasta_file = sys.argv[1]
length_threshold = int(sys.argv[2])

filter_fasta_by_length(fasta_file, length_threshold)
