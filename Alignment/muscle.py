import os
import itertools
import subprocess
import sys
from Bio import SeqIO

def write_pair_to_fasta(seq1_name, seq1_data, seq2_name, seq2_data, temp_fasta_path):
    """Writes a pair of sequences to a FASTA file."""
    with open(temp_fasta_path, 'w') as fasta_file:
        fasta_file.write(f">{seq1_name}\n{seq1_data}\n")
        fasta_file.write(f">{seq2_name}\n{seq2_data}\n")

def create_slurm_script(batch_number, batch_pairs, output_dir, job_name):
    """Generates a SLURM job script for a batch of sequence pairs."""
    batch_dir = os.path.join(output_dir, f"batch_{batch_number}")
    os.makedirs(batch_dir, exist_ok=True)

    slurm_script = os.path.join(batch_dir, f"job_{batch_number}.sh")
    with open(slurm_script, 'w') as slurm_file:
        slurm_file.write("#!/bin/bash\n")
        slurm_file.write(f"#SBATCH --job-name={job_name}_{batch_number}\n")
        slurm_file.write("#SBATCH --time=04:00:00\n")
        slurm_file.write("#SBATCH --ntasks=1\n")
        slurm_file.write("#SBATCH --cpus-per-task=1\n")
        slurm_file.write("#SBATCH --mem=3G\n")
        slurm_file.write("module load MUSCLE/3.8.1551-GCC-9.3.0\n") 
        slurm_file.write("cd $SLURM_SUBMIT_DIR\n")
        
        for pair in batch_pairs:
            seq1, seq2 = pair  # Unpack the pair into two sequences (not four)
            seq1_name, seq1_data = seq1.id, str(seq1.seq)
            seq2_name, seq2_data = seq2.id, str(seq2.seq)
            
            temp_fasta_path = os.path.join(batch_dir, f"{seq1_name}_{seq2_name}.fasta")
            
            # Write the FASTA pair file
            write_pair_to_fasta(seq1_name, seq1_data, seq2_name, seq2_data, temp_fasta_path)
            
            # Add command to run muscle on the generated FASTA
            slurm_file.write(f"muscle -in {temp_fasta_path} -out {temp_fasta_path.replace('.fasta', '.aln')}\n")
            
            # Clean up FASTA file after processing
            slurm_file.write(f"rm {temp_fasta_path}\n")
        
        # Optional: Add more commands for post-processing if necessary.

    # Submit the SLURM script
    subprocess.run(["sbatch", slurm_script])

def generate_jobs(fasta_file, output_dir, batch_size):
    """Generates pairs from the input FASTA and creates SLURM job scripts."""
    # Read the sequences from the input FASTA
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    
    # Generate all possible pairs (combinations of 2)
    sequence_pairs = list(itertools.combinations(sequences, 2))
    
    # Split the pairs into batches
    batches = [sequence_pairs[i:i + batch_size] for i in range(0, len(sequence_pairs), batch_size)]
    
    # Generate job scripts for each batch
    for batch_number, batch_pairs in enumerate(batches, start=1):
        create_slurm_script(batch_number, batch_pairs, output_dir, "muscle_job")

def main():
    """Main function to run the job generation."""
    # Parse command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python3 ./muscle.py <input_fasta> <output_dir> <batch_size>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]  # Input FASTA file
    output_dir = sys.argv[2]  # Directory to save batch files and job scripts
    batch_size = int(sys.argv[3])  # Number of pairs per batch
    
    # Generate jobs
    generate_jobs(fasta_file, output_dir, batch_size)

if __name__ == "__main__":
    main()
