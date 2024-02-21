This repository describes how the analyses for the Lepidoptera pherome paper were performed. 

# blast

The Blast folder contains all the scripts used to:
1) Retrieve all the CDSs from *Heliconius melpomene* v2.5 reference genome (lepbase).
2) Retrieve FAR or FAD sequences from uniprot and genbank.
3) Create a blast database using the FAR or FAD sequences.
4) Blast(x) each *Heliconius melpomene* CDS sequence against the FAR of the FAD database.
5) Extract the best blastx hits
6) Check the location in the *Heliconius melpomene* reference genome of best hit.
7) Plot the result in R using the package chromoMap

To run the script just run:
```
./master.sh FAR
./master.sh FAD
```

This script requires [Biopython](http://biopython.org/), [edirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/), and [chromoMap](https://cran.r-project.org/web/packages/chromoMap/vignettes/chromoMap.html). 
