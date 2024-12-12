#!/usr/bin/python
from Bio import SeqIO
import sys
import os

##################################################################
# Command that look at the size of each sequence in a fasta file #
##################################################################
#usage: python -W ignore *fasta

input_seq_iterator = SeqIO.parse(open(sys.argv[1], "r"), "fasta")

for record in input_seq_iterator:
	print(record.id + "\t" + str(1) + "\t"  + str(len(record.seq)))
