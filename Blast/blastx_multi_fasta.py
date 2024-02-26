#!/usr/bin/env python

import sys
import os
import subprocess
from Bio import SeqIO

#############################
# Function that runs blastx #
#############################
def run_blastx(query_sequence, db_path):
	blastx_cmd = [
		'blastx',
		'-query', '-',
		'-db', db_path,
		'-outfmt', '7'
	]

	try:
		input_bytes = query_sequence   # Encode using UTF-8
		result = subprocess.run(blastx_cmd, input=input_bytes, text=True, capture_output=True, check=True)
		return result.stdout
	except subprocess.CalledProcessError as e:
		print(f"Error occurred during BLASTP search: {e}")
		return None

############################################################
# Function that get the blastx results in the right format #
############################################################
def main(input_file, array_task_id, db_path):
	# Extract directory path and filename from input file
	input_dir, input_filename = os.path.split(input_file)
	# Construct output file path using input directory and array task ID
	output_file_name = os.path.join(input_dir, f"blast_hits_{array_task_id}.txt")
	
	with open(output_file_name, "w") as blast_output_file:
		blast_output_file.write("query\tsubject\tidentity\talignmentlength\tmismatches\tgapopens\tquerystart\tqueryend\tsubjectstart\tsubjectend\tevalue\tbitscore\n")
		for record in SeqIO.parse(input_file, "fasta"):
			# Extract sequence identifier from the header
			query_name = record.id.split()[0]  # Take the part before the first whitespace
			query_sequence = str(record.seq)
			print(f"Running BLASTP for sequence: {query_name}")
			blastx_result = run_blastx(query_sequence,db_path)
			if "# 0 hits found" not in blastx_result:
				for line in blastx_result.splitlines():
					if not line.startswith("#"):
						# Split the line by tab and omit the second column (Query_1)
						columns = line.split("\t")
						modified_line = "\t".join(columns[1:])
						if(float(columns[-1]) > 100):
							blast_output_file.write(f"{query_name}\t{modified_line}\n")


if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("Usage: python blastx_multi_fasta.py input.fasta array_task_id db_path")
		sys.exit(1)
	input_file = sys.argv[1]
	array_task_id = sys.argv[2]
	db_path = sys.argv[3]
	main(input_file, array_task_id, db_path)
