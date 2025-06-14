#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def parse_gff(gff_file):
    # Utiliser (seqid, gene_start, gene_end) comme clé unique pour chaque gene
    cds_dict = defaultdict(list)
    gene_ids = {}
    current_gene = None
    with open(gff_file) as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 8:
                continue
            if len(parts) == 8:
                parts.append(".")
            seqid, source, feature, start, end, score, strand, phase, attributes = parts
            
            if feature == "gene":
                # Construire clé unique gene_key = (seqid, start, end)
                gene_key = (seqid, int(start), int(end))
                # Extraire un identifiant unique depuis attributes (tu peux adapter selon ton format)
                gid = attributes.split("sequence")[-1].strip().replace(";", "").replace(" ", "_")
                gene_ids[gene_key] = gid
                current_gene = gene_key
            
            elif feature == "cds" and current_gene is not None:
                # Ajouter le CDS à ce gène (clé current_gene)
                cds_dict[current_gene].append((int(start), int(end), strand))
    return cds_dict, gene_ids

def trim_exon(seq):
    length = len(seq)
    trim_len = length % 3
    if trim_len != 0:
        seq = seq[:-trim_len]
    return seq

def extract_cds(cds_dict, gene_ids, genome_fasta, output_cds, output_prot):
    genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
    with open(output_cds, "w") as out_cds, open(output_prot, "w") as out_prot:
        for gene_key in cds_dict:
            seqid = gene_key[0]
            if seqid not in genome:
                continue
            scaffold_seq = genome[seqid].seq
            strand = cds_dict[gene_key][0][2]
            exons = sorted(cds_dict[gene_key], key=lambda x: x[0], reverse=(strand == "-"))
            full_seq = ""
            for start, end, _ in exons:
                exon_seq = scaffold_seq[start - 1:end]
                exon_seq = trim_exon(exon_seq)
                full_seq += exon_seq
            if strand == "-":
                full_seq = full_seq.reverse_complement()
            gene_id = gene_ids.get(gene_key, f"{seqid}_{gene_key[1]}_{gene_key[2]}")
            # Write CDS
            out_cds.write(f">{gene_id}\n{full_seq}\n")
            # Translate and write protein
            protein_seq = full_seq.translate(to_stop=True)
            out_prot.write(f">{gene_id}\n{protein_seq}\n")

def main():
    parser = argparse.ArgumentParser(description="Extract CDS from Exonerate GFF with exon trimming")
    parser.add_argument("-f", "--gff", required=True, help="Exonerate GFF file")
    parser.add_argument("-g", "--genome", required=True, help="Genome FASTA file")
    parser.add_argument("-o", "--output", default="cds.fasta", help="Output CDS FASTA file prefix")
    args = parser.parse_args()

    cds_dict, gene_ids = parse_gff(args.gff)
    extract_cds(cds_dict, gene_ids, args.genome, args.output, args.output.replace(".fasta", "_prot.fasta"))

if __name__ == "__main__":
    main()
