#!/usr/bin/env python

import sys  # Importing the sys module to access command-line arguments
import os  # Importing the os module to perform operating system related tasks
import subprocess  # Importing the subprocess module to run external commands
from Bio import SeqIO  # Importing SeqIO module from Biopython for sequence parsing

# Function to run BLASTX with given query sequence and database path
def run_blastx(query_sequence, db_path):
	# BLASTX command and its arguments
	blastx_cmd = [
		'blastx',  # BLASTX program name
		'-query', '-',  # Use standard input as query sequence
		'-db', db_path,  # Path to the BLAST database
		'-outfmt', '7'  # Output format
	]

	try:
		# Encode query sequence using UTF-8
		input_bytes = query_sequence
		# Run the BLASTX command, capturing the output
		result = subprocess.run(blastx_cmd, input=input_bytes, text=True, capture_output=True, check=True)
		# Return the stdout of the BLASTX command
		return result.stdout
	except subprocess.CalledProcessError as e:
		# Handle error if BLASTX command returns non-zero exit status
		print(f"Error occurred during BLASTP search: {e}")
		return None

# Main function to perform BLASTX search on multiple sequences
def main(input_file, array_task_id, db_path):
	# Extract directory path and filename from input file
	input_dir, input_filename = os.path.split(input_file)
	# Construct output file path using input directory and array task ID
	output_file_name = os.path.join(input_dir, f"blast_hits_{array_task_id}.txt")
	
	# Open output file in write mode
	with open(output_file_name, "w") as blast_output_file:
		# Write header to output file
		blast_output_file.write("query\tsubject\tidentity\talignmentlength\tmismatches\tgapopens\tquerystart\tqueryend\tsubjectstart\tsubjectend\tevalue\tbitscore\n")
		# Iterate over each sequence in the input FASTA file
		for record in SeqIO.parse(input_file, "fasta"):
			# Extract sequence identifier from the header
			query_name = record.id.split()[0]  # Take the part before the first whitespace
			# Convert sequence to string
			query_sequence = str(record.seq)
			# Print message indicating BLASTP is running for current sequence
			print(f"Running BLASTP for sequence: {query_name}")
			# Run BLASTX for current sequence and get the result
			blastx_result = run_blastx(query_sequence, db_path)
			# Check if any hits were found in BLASTX result
			if "# 0 hits found" not in blastx_result:
				# Iterate over each line in BLASTX result
				for line in blastx_result.splitlines():
					# Check if the line is not a comment line
					if not line.startswith("#"):
						# Split the line by tab and omit the second column (Query_1)
						columns = line.split("\t")
						modified_line = "\t".join(columns[1:])  # Concatenate all columns except the first one
						# Check if bitscore is greater than 100
						if float(columns[-1]) > 100:  # Convert the last column (bitscore) to float and check if it's greater than 100
							# Write query name and modified line to output file
							blast_output_file.write(f"{query_name}\t{modified_line}\n")

# Check if the script is being run as the main program
if __name__ == "__main__":
	# Check if the correct number of command-line arguments is provided
	if len(sys.argv) != 4:
		# Print usage message and exit with status 1 if incorrect number of arguments
		print("Usage: python blastx_multi_fasta.py input.fasta array_task_id db_path")
		sys.exit(1)
	# Get input file path, array task ID, and database path from command-line arguments
	input_file = sys.argv[1]
	array_task_id = sys.argv[2]
	db_path = sys.argv[3]
	# Call the main function with provided arguments
	main(input_file, array_task_id, db_path)
