args <- commandArgs(trailingOnly = TRUE)
alignment <- args[1]
output <- args[2]

# Load libraries
library(msaR)
library(htmlwidgets)
library(Biostrings)

alignm <- Biostrings::readAAStringSet(alignment)
print(alignm)

# Create the msaR object
msa_output <- msaR(alignm)
print(msa_output)

# Save the msaR output as an HTML file
saveWidget(msa_output, output, selfcontained = TRUE)
