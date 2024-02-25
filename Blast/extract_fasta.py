import sys
from Bio import SeqIO

def extract_sequences(id_file, fasta_file, output_file):
	# Read IDs from the file
	with open(id_file, 'r') as f:
		ids = {line.strip() for line in f}

	# Extract sequences from the FASTA file
	with open(output_file, 'w') as out_f:
		for record in SeqIO.parse(fasta_file, 'fasta'):
			if record.id in ids:
				SeqIO.write(record, out_f, 'fasta')

if __name__ == "__main__":
	if len(sys.argv) != 4:
		print("Usage: python extract_fasta.py <id_file> <fasta_file> <output_file>")
		sys.exit(1)
	
	id_file = sys.argv[1]
	fasta_file = sys.argv[2]
	output_file = sys.argv[3]

	extract_sequences(id_file, fasta_file, output_file)
