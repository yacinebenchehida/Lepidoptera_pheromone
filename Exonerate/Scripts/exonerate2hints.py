# save as exonerate2hints.py
import sys

def convert_exonerate_to_hints(infile, outfile):
    with open(infile) as fin, open(outfile, "w") as fout:
        for line in fin:
            if line.startswith("#"): continue
            parts = line.strip().split("\t")
            if len(parts) < 9: continue
            source = parts[1]
            feature = parts[2]
            if feature not in ["exon", "cds"]: continue  # use exon/CDS only
            start = parts[3]
            end = parts[4]
            strand = parts[6]
            fout.write(f"{parts[0]}\t{source}\t{feature}\t{start}\t{end}\t.\t{strand}\t.\tsrc=E;pri=4;\n")

if __name__ == "__main__":
    convert_exonerate_to_hints(sys.argv[1], sys.argv[2])
