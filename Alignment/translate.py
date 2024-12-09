from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

def translate_fasta(handle):
    """Translate nucleotide sequences to proteins and print to stdout in FASTA format."""
    for record in SeqIO.parse(handle, "fasta"):
        translated_seq = record.seq.translate(to_stop=True)
        translated_record = SeqRecord(translated_seq, id=record.id, description="")  # No description
        SeqIO.write(translated_record, sys.stdout, "fasta")

if __name__ == "__main__":
    if len(sys.argv) == 2:
        input_fasta = sys.argv[1]
        with open(input_fasta, "r") as handle:
            translate_fasta(handle)
    elif not sys.stdin.isatty():  # Check if input is provided via pipe
        translate_fasta(sys.stdin)
    else:
        print("Usage: python3 translate_fasta.py <input_fasta> OR provide input via pipe", file=sys.stderr)
        sys.exit(1)
