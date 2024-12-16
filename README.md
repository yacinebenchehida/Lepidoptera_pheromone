This repository describes how the analyses for the Lepidoptera pheromone paper were performed. 

# blast

The blast folder includes scripts to:

1) Retrieve all coding sequences (CDS) for several Lepidoptera species from LepBase and ENSEMBL.
2) Retrieve FAR and FAD sequences from UniProt and GenBank.
3) Create a BLAST database using the FAR or FAD sequences.
4) Perform BLASTx of the CDS sequences of each species against the FAR or FAD database.
5) Extract the best BLASTx hits.
6) Determine the genomic location of the best hits in each species' reference genome.
7) Visualize the results as ideograms for each species.

# Alignment

The blast folder includes scripts to:

1) Identify paralogs based on sequence identity. 
2) Construct phylogenies of all FAR/FAD genes across species.

# How to run the script

Run the main script as follows:
```
./master.sh FAR FAD 
```

This script requires [Biopython](http://biopython.org/) and [edirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (need to update dependencies)
