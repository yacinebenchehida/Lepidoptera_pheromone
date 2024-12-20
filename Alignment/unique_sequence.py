import sys
from Bio import SeqIO

def remove_duplicate_sequences(input_fasta, output_fasta):
    seen_sequences = set()  # To track unique sequences
    unique_records = []     # To store the final records
    
    # Open the input FASTA file and read the sequences
    with open(input_fasta, "r") as input_file:
        for record in SeqIO.parse(input_file, "fasta"):
            # If this sequence hasn't been seen before, add it to the unique list
            if str(record.seq) not in seen_sequences:
                seen_sequences.add(str(record.seq))
                unique_records.append(record)
    
    # Write the unique sequences back to the output file
    with open(output_fasta, "w") as output_file:
        SeqIO.write(unique_records, output_file, "fasta")
    print(f"Filtered FASTA file saved to {output_fasta}")

# Ensure that the script is run with input and output arguments
if len(sys.argv) != 3:
    print("Usage: python script.py <input_fasta_file> <output_fasta_file>")
    sys.exit(1)

# Read the input and output file paths from command-line arguments
input_fasta = sys.argv[1]
output_fasta = sys.argv[2]

# Run the function
remove_duplicate_sequences(input_fasta, output_fasta)
