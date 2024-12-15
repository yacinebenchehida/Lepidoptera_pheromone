This repository describes how the analyses for the Lepidoptera pheromone paper were performed. 

# blast

The Blast folder contains all the scripts used to:
1) Retrieve all the CDSs from several lepidopteran species from lepbase and ENSEMBL.
2) Retrieve FAR or FAD sequences from uniprot and genbank.
3) Create a blast database using the FAR or FAD sequences.
4) Blast(x) the CDS sequence of each species against the FAR of the FAD database.
5) Extract the best blastx hits
6) Check the location in each species reference genome of best hits.
7) Plot the result as an ideogram for each species.
8) Make a phylogeny of all FAR/FAD genes for all species.

To run the script just run:
```
./master.sh FAR FAD
```

This script requires [Biopython](http://biopython.org/) and [edirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/) (need to update dependencies)
