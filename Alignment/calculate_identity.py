import os
import sys
from Bio import AlignIO

def calculate_identity_and_stretch(alignment_file):
    """Calculate the percentage identity and gap statistics for the alignment."""
    # Read the alignment using BioPython
    alignment = AlignIO.read(alignment_file, "fasta")
    
    # Ensure there are exactly two sequences in the alignment
    seq1, seq2 = alignment[0], alignment[1]
    
    total_residues = 0
    matching_residues = 0
    seqs = [seq1, seq2]

    # Calculate percentage identity
    for i in range(len(seq1)):
        if seq1[i] != "-" and seq2[i] != "-":  # Exclude gaps
            total_residues += 1
            if seq1[i] == seq2[i]:  # Check for identical residues
                matching_residues += 1
    
    # Compute the percentage identity
    identity = (matching_residues / total_residues) * 100 if total_residues > 0 else 0

    # Identify the sequence with the most gaps
    max_gaps_seq = max(seqs, key=lambda seq: seq.seq.count('-'))
    total_length = len(max_gaps_seq.seq)
    num_gaps = max_gaps_seq.seq.count('-')
    
    # Compute percentage of gaps in the gappiest sequence
    gap_percentage = (num_gaps / total_length) * 100 if total_length > 0 else 0

    # Compute the longest stretch without gaps and count stretches
    max_stretch = 0
    current_stretch = 0
    stretch_count = 0
    in_stretch = False  # Track whether we are in a stretch
    
    for base in max_gaps_seq.seq:
        if base != "-":
            current_stretch += 1
            max_stretch = max(max_stretch, current_stretch)
            if not in_stretch:
                stretch_count += 1
                in_stretch = True
        else:
            current_stretch = 0
            in_stretch = False
    
    # Compute max_stretch ratio
    non_gap_count = total_length - num_gaps
    max_stretch_ratio = (max_stretch / non_gap_count) if non_gap_count > 0 else 0

    # Return results: identity, sequence names, max stretch, gap percentage, stretch count, and max_stretch_ratio
    return identity, seq1.id, seq2.id, max_stretch, gap_percentage, stretch_count, max_stretch_ratio

def process_batch(batch_folder):
    """Process all alignment files in a batch folder and calculate identity and gap statistics."""
    # Iterate through all files in the batch folder
    for file_name in os.listdir(batch_folder):
        if file_name.endswith(".aln"):  # Process only .aln files
            alignment_file = os.path.join(batch_folder, file_name)
            identity, seq1_name, seq2_name, max_stretch, gap_percentage, stretch_count, max_stretch_ratio = calculate_identity_and_stretch(alignment_file)
            
            # Output file with the same name as the alignment but with .txt extension
            identity_output = alignment_file.replace(".aln", ".txt")
            
            # Write results to the output text file
            with open(identity_output, "w") as f:
                # Write results in the format:
                # seq1 seq2 percentage_identity max_stretch max_stretch_ratio gap_percentage stretch_count
                f.write(f"{seq1_name} {seq2_name} {identity:.2f} {max_stretch} {max_stretch_ratio:.4f} {gap_percentage:.2f} {stretch_count}\n")

def main():
    """Main function to process all batch folders."""
    # Check input arguments
    if len(sys.argv) != 2:
        print("Usage: python3 calculate_identity.py <batch_directory>")
        sys.exit(1)
    
    # Directory containing batch subdirectories
    batch_dir = sys.argv[1]
    
    # Verify that the directory exists
    if not os.path.isdir(batch_dir):
        print(f"Error: The directory '{batch_dir}' does not exist.")
        sys.exit(1)

    # Iterate through all batch subdirectories
    for batch_folder in os.listdir(batch_dir):
        batch_folder_path = os.path.join(batch_dir, batch_folder)
        
        if os.path.isdir(batch_folder_path):
            print(f"Processing batch folder: {batch_folder_path}")
            process_batch(batch_folder_path)

if __name__ == "__main__":
    main()

