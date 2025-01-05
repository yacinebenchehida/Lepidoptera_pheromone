import re
import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq

def extract_unique_sequences(input_fasta, output_fasta):
    # Dictionary to store species as keys and a list of sequences (with their headers) as values
    species_sequences = {}

    # Parse the FASTA file
    for record in SeqIO.parse(input_fasta, "fasta"):
        header = record.id
        sequence = str(record.seq)
        
        # Extract species using the updated regex (genus_species format)
        match = re.search(r"([A-Za-z]+_[A-Za-z]+)$", header)
        if not match:
            raise ValueError(f"Header format incorrect for: {header}")
        
        species = match.group(1)
        
        # Initialize species entry if not present
        if species not in species_sequences:
            species_sequences[species] = {}

        # Store the sequence with its header, ensuring unique sequences per species
        if sequence not in species_sequences[species]:
            species_sequences[species][sequence] = header

    # Prepare records for output
    unique_records = []
    for species, seq_map in species_sequences.items():
        for sequence, header in seq_map.items():
            # For each species, we keep the original sequence with its header
            unique_records.append(SeqIO.SeqRecord(seq=Seq(sequence), id=header, description=""))

    # Write the unique sequences to the output file
    SeqIO.write(unique_records, output_fasta, "fasta")

if __name__ == "__main__":
    # Ensure proper arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    extract_unique_sequences(input_fasta, output_fasta)
