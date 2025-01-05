import sys
from Bio import SeqIO
import argparse

def filter_sequences(input_handle, output_file, kept_samples_file, ratio, max_ratio):
    # Read all sequences
    sequences = list(SeqIO.parse(input_handle, "fasta"))

    # Calculate the mean length of all sequences
    mean_length = sum(len(record.seq) for record in sequences) / len(sequences)

    # Calculate the length thresholds based on the ratios
    min_threshold = mean_length * ratio
    max_threshold = mean_length * max_ratio
    print(f"Mean sequence length: {mean_length:.2f}")
    print(f"Minimum length threshold ({ratio}): {min_threshold:.2f}")
    print(f"Maximum length threshold ({max_ratio}): {max_threshold:.2f}")

    # Filter sequences
    filtered_sequences = [
        record for record in sequences 
        if len(record.seq) >= min_threshold and len(record.seq) <= max_threshold
    ]

    print(f"Retained {len(filtered_sequences)} out of {len(sequences)} sequences.")

    # Write the filtered sequences to the output file
    with open(output_file, "w") as out_handle:
        SeqIO.write(filtered_sequences, out_handle, "fasta")

    # Write the list of kept sample IDs to a single-line file
    kept_samples = " ".join(record.id for record in filtered_sequences)
    with open(kept_samples_file, "w") as sample_handle:
        sample_handle.write(kept_samples + "\n")

def main():
    parser = argparse.ArgumentParser(description="Filter sequences based on length from a FASTA file.")
    parser.add_argument(
        "-i", "--input", 
        type=str, 
        default=None,
        help="Input FASTA file. Use '-' to read from stdin."
    )
    parser.add_argument(
        "-o", "--output", 
        type=str, 
        required=True,
        help="Output FASTA file with filtered sequences."
    )
    parser.add_argument(
        "-s", "--samples", 
        type=str, 
        required=True,
        help="Output file to write the list of kept sample IDs (space-separated in one line)."
    )
    parser.add_argument(
        "-r", "--ratio", 
        type=float, 
        default=0.8,
        help="Minimum sequence length as a fraction of the mean length (default: 0.8)."
    )
    parser.add_argument(
        "-m", "--max-ratio", 
        type=float, 
        default=1.2,
        help="Maximum sequence length as a fraction of the mean length (default: 1.2)."
    )

    args = parser.parse_args()

    # Use stdin if input file is not provided or explicitly set to "-"
    if args.input == "-" or args.input is None:
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