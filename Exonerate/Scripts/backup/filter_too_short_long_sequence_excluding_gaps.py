import sys
from Bio import SeqIO
import argparse

def ungapped_length(seq):
    return len(str(seq).replace("-", "").replace(".", ""))

def filter_sequences(input_handle, output_file, kept_samples_file, ratio, max_ratio):
    # Read all sequences
    sequences = list(SeqIO.parse(input_handle, "fasta"))

    # Calculate ungapped lengths
    ungapped_lengths = [ungapped_length(record.seq) for record in sequences]

    # Calculate mean ungapped length
    mean_length = sum(ungapped_lengths) / len(ungapped_lengths)

    # Calculate thresholds
    min_threshold = mean_length * ratio
    max_threshold = mean_length * max_ratio

    print(f"Mean ungapped sequence length: {mean_length:.2f}")
    print(f"Minimum ungapped length threshold ({ratio}): {min_threshold:.2f}")
    print(f"Maximum ungapped length threshold ({max_ratio}): {max_threshold:.2f}")

    # Filter sequences based on ungapped length
    filtered_sequences = [
        record for record, ungapped_len in zip(sequences, ungapped_lengths)
        if min_threshold <= ungapped_len <= max_threshold
    ]

    print(f"Retained {len(filtered_sequences)} out of {len(sequences)} sequences.")

    # Write filtered sequences
    with open(output_file, "w") as out_handle:
        SeqIO.write(filtered_sequences, out_handle, "fasta")

    # Write retained sample IDs
    kept_samples = " ".join(record.id for record in filtered_sequences)
    with open(kept_samples_file, "w") as sample_handle:
        sample_handle.write(kept_samples + "\n")

def main():
    parser = argparse.ArgumentParser(description="Filter sequences based on ungapped length from a FASTA file.")
    parser.add_argument("-i", "--input", type=str, default=None,
                        help="Input FASTA file. Use '-' to read from stdin.")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Output FASTA file with filtered sequences.")
    parser.add_argument("-s", "--samples", type=str, required=True,
                        help="Output file to write the list of kept sample IDs (space-separated).")
    parser.add_argument("-r", "--ratio", type=float, default=0.8,
                        help="Minimum ungapped length as a fraction of the mean (default: 0.8).")
    parser.add_argument("-m", "--max-ratio", type=float, default=1.2,
                        help="Maximum ungapped length as a fraction of the mean (default: 1.2).")
    args = parser.parse_args()

    # Input from stdin or file
    if args.input in [None, "-"]:
        input_handle = sys.stdin
    else:
        input_handle = open(args.input, "r")

    try:
        filter_sequences(input_handle, args.output, args.samples, args.ratio, args.max_ratio)
    finally:
        if args.input not in [None, "-"]:
            input_handle.close()

if __name__ == "__main__":
    main()