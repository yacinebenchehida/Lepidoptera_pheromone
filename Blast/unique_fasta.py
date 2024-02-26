#!/usr/bin/env python

import sys
from collections import defaultdict
from Bio import SeqIO
import os

# Function to print unique sequences to the screen
def print_unique_sequences(input_file, output_folder):
	# Initialize dictionaries to store sequences and their headers
	sequences = {}
	headers = defaultdict(list)

	# Open the FASTA file and iterate over each record
	with open(input_file, "r") as f:
		for record in SeqIO.parse(f, "fasta"):
			# Get the sequence and its header
			sequence = str(record.seq)
			header = record.description

			# Check if the sequence is already present
			if sequence in sequences:
				# If so, append the header to the list of headers
				headers[sequence].append(header)
			else:
				# Otherwise, store the sequence and its header
				sequences[sequence] = header

	# Create the output folder if it does not exist
	os.makedirs(output_folder, exist_ok=True)
	non_unique_file_path = os.path.join(output_folder, "non_unique.txt")

	# Print unique sequences to the screen and store non-unique headers
	with open(non_unique_file_path, "w") as non_unique_file:
		for i, (sequence, header) in enumerate(sequences.items(), start=1):
			PREFIX = sys.argv[2]
			print(">" + PREFIX + "_" + str(i))
			print(sequence)
			if sequence in headers:
				non_unique_file.write("-------" + PREFIX + "_{}".format(i))
				non_unique_file.write("-------\n")
				non_unique_file.write(header + "\n")
				for non_unique_header in headers[sequence]:
					non_unique_file.write(non_unique_header + "\n")
		non_unique_file.write("--------------\n")

# Main function
def main():
	# Check if the number of command-line arguments is correct
	if len(sys.argv) != 3:
		# If not, print usage information and exit with error code 1
		print("Usage: python ./unique_fasta.py input.fasta prefix")
		sys.exit(1)

	# Extract the input file path and prefix from the command-line arguments
	input_file = sys.argv[1]
	output_folder = os.path.dirname(input_file)

	# Print unique sequences from the input FASTA file
	print_unique_sequences(input_file, output_folder)

# Entry point of the script
if __name__ == "__main__":
	# Call the main function when the script is executed
	main()
