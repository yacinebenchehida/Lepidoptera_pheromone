import os
import re
import argparse
from Bio import SeqIO

def extract_genus_species(header):
    """
    Extract Genus_Species from a fasta header.
    """
    match = re.search(r"^[^_]+_([^_]+_[^_]+)$", header)
    return match.group(1) if match else None

def process_fasta(file_path):
    """
    Process a fasta file to check for at least 3 distinct Genus_Species.
    Returns a list of headers if valid, else None.
    """
    genus_species_set = set()
    headers = []

    # Read the fasta file
    for record in SeqIO.parse(file_path, "fasta"):
        header = record.description
        genus_species = extract_genus_species(header)
        if genus_species:
            genus_species_set.add(genus_species)
            headers.append(header)

    # Check if there are at least 3 unique Genus_Species
    if len(genus_species_set) >= 3:
        return headers
    return None

def main(folder_path, pheromone_name):
    """
    Main function to identify and process valid fasta files.
    """
    # Define output folder
    output_folder = os.path.join(folder_path, f"{pheromone_name}_manual_curation")
    os.makedirs(output_folder, exist_ok=True)
    output_file = os.path.join(output_folder, f"Merged_manual_curation_{pheromone_name}.txt")

    merged_headers = []  # Store headers from "good" fasta files

    # Loop through numbered folders starting with pheromone_name
    folder_number = 1
    while True:
        subfolder = os.path.join(folder_path, f"{pheromone_name}_{folder_number}")

        if not os.path.isdir(subfolder):
            # Stop when a folder doesn't exist
            break

        # Look for files ending with "nt_trimmed.aln" in the subfolder
        for filename in os.listdir(subfolder):
            if filename.endswith("nt_trimmed_unique_seq_sp.aln"):
                fasta_path = os.path.join(subfolder, filename)
                print(f"Processing: {fasta_path}")
                
                # Process the fasta file
                headers = process_fasta(fasta_path)
                if headers:
                    print(f"File {filename} is valid. Adding headers.")
                    merged_headers.extend(headers)
                else:
                    print(f"File {filename} is invalid (less than 3 unique Genus_Species).")

        folder_number += 1

    # Write all headers to the output file
    if merged_headers:
        with open(output_file, "w") as out_f:
            out_f.write(" ".join(merged_headers))
        print(f"Headers merged and saved to: {output_file}")
    else:
        print("No valid fasta files found.")

if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process FASTA files for manual curation.")
    parser.add_argument("folder_path", help="Path to the folder containing pheromone subfolders")
    parser.add_argument("pheromone_name", help="Name of the pheromone gene")
    args = parser.parse_args()

    # Run the main function
    main(args.folder_path, args.pheromone_name)

