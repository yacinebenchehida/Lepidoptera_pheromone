import sys
import os
import subprocess
from Bio import SeqIO

def run_blastp(query_sequence):
	blastp_cmd = [
		'blastp',
		'-query', '-',
		'-db', '/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results/FAR_db',
		'-outfmt', '7'
	]

	try:
		input_bytes = query_sequence   # Encode using UTF-8
		result = subprocess.run(blastp_cmd, input=input_bytes, text=True, capture_output=True, check=True)
		return result.stdout
	except subprocess.CalledProcessError as e:
		print(f"Error occurred during BLASTP search: {e}")
		return None

def main(input_file, array_task_id):
	# Extract directory path and filename from input file
	input_dir, input_filename = os.path.split(input_file)
	# Construct output file path using input directory and array task ID
	output_file_name = os.path.join(input_dir, f"blast_FAR_hits_{array_task_id}.txt")
	
	with open(output_file_name, "w") as blast_output_file:
		blast_output_file.write("query\tsubject\tidentity\talignmentlength\tmismatches\tgapopens\tquerystart\tqueryend\tsubjectstart\tsubjectend\tevalue\tbitscore\n")
		for record in SeqIO.parse(input_file, "fasta"):
			# Extract sequence identifier from the header
			query_name = record.id.split()[0]  # Take the part before the first whitespace
			query_sequence = str(record.seq)
			print(f"Running BLASTP for sequence: {query_name}")
			blastp_result = run_blastp(query_sequence)
			if "# 0 hits found" not in blastp_result:
				for line in blastp_result.splitlines():
					if not line.startswith("#"):
						# Split the line by tab and omit the second column (Query_1)
						columns = line.split("\t")
						modified_line = "\t".join(columns[1:])
						if(float(columns[-1]) > 50):
							blast_output_file.write(f"{query_name}\t{modified_line}\n")

if __name__ == "__main__":
	if len(sys.argv) != 3:
		print("Usage: python blastp_multi_fasta.py input.fasta array_task_id")
		sys.exit(1)
	input_file = sys.argv[1]
	array_task_id = sys.argv[2]
	main(input_file, array_task_id)
