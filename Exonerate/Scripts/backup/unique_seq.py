import argparse
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Print unique sequences from a FASTA file.")
    parser.add_argument("input_fasta", help="Input FASTA file")
    args = parser.parse_args()

    seen_seqs = set()
    unique_seqs = []

    for record in SeqIO.parse(args.input_fasta, "fasta"):
        seq_str = str(record.seq)
        if seq_str not in seen_seqs:
            seen_seqs.add(seq_str)
            unique_seqs.append(seq_str)

    for i, seq in enumerate(unique_seqs, 1):
        print(f">FAD_{i}")
        print(seq)

if __name__ == "__main__":
    main()
