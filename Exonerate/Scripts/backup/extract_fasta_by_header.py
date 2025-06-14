from Bio import SeqIO
import sys

def get_sequences_by_headers(fasta_file, headers):
    """
    Given a fasta file and a list of headers, this function returns the matching sequences.
    
    Parameters:
    fasta_file (str): Path to the input FASTA file.
    headers (list): List of FASTA headers to extract sequences for.
    
    Returns:
    None: It prints the headers and corresponding sequences that match the provided headers.
    """
    
    try:
        # Read the FASTA file
        with open(fasta_file, "r") as handle:
            # Iterate through the sequences in the FASTA file
            for record in SeqIO.parse(handle, "fasta"):
                # Check if the header matches any of the provided headers
                if record.id in headers:
                    # Print the header and corresponding sequence
                    print(f">{record.id}")
                    print(record.seq)
    except FileNotFoundError:
        print(f"Error: The file {fasta_file} was not found.")
    except Exception as e:
        print(f"An error occurred: {str(e)}")


def read_headers_from_file(header_file):
    """
    Read headers from a file where headers are space-separated in a single line.
    """
    try:
        with open(header_file, 'r') as file:
            headers = file.read().strip().split()
        return headers
    except Exception as e:
        print(f"Error reading headers from file: {e}")
        return []


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python script.py <fasta_file> --file <headers_file> OR --list <header1/header2...>")
        sys.exit(1)

    fasta_file = sys.argv[1]
    input_type = sys.argv[2]
    headers_input = sys.argv[3:]

    if input_type == '--file':
        headers = read_headers_from_file(headers_input[0])
    elif input_type == '--list':
        headers = headers_input
    else:
        print("Invalid input type. Use --file for header file or --list for list of headers.")
        sys.exit(1)

    # Call the function to retrieve and print the sequences
    get_sequences_by_headers(fasta_file, headers)