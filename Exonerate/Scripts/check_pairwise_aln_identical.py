from Bio import SeqIO
from itertools import combinations

def load_alignment(fasta_file):
    return list(SeqIO.parse(fasta_file, "fasta"))

def aligned_identity(seq1, seq2):
    aligned1 = str(seq1.seq)
    aligned2 = str(seq2.seq)
    
    if len(aligned1) != len(aligned2):
        return False

    # Compare only aligned (non-gap in both) positions
    for a, b in zip(aligned1, aligned2):
        if a != "-" and b != "-":
            if a != b:
                return False
    return True

def group_identical(records):
    groups = []
    visited = set()

    for i, rec1 in enumerate(records):
        if i in visited:
            continue
        group = [rec1]
        visited.add(i)
        for j in range(i + 1, len(records)):
            if j in visited:
                continue
            rec2 = records[j]
            if aligned_identity(rec1, rec2):
                group.append(rec2)
                visited.add(j)
        groups.append(group)
    return groups

def print_longest_per_group(groups):
    for group in groups:
        # Longest full aligned sequence (with gaps)
        longest = max(group, key=lambda r: len(str(r.seq)))
        print(f">{longest.id}\n{longest.seq}")

# Example usage
if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        print("Usage: python script.py alignment.fasta")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    records = load_alignment(fasta_file)
    groups = group_identical(records)
    print_longest_per_group(groups)
