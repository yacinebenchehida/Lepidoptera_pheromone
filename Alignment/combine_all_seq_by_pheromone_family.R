library(Biostrings)
library(tibble)
library(dplyr)

# Function to pad sequences with gaps
pad_with_gaps <- function(seq, length_needed) {
  current_length <- nchar(as.character(seq))
  if (current_length < length_needed) {
    paste0(seq, strrep("-", length_needed - current_length))
  } else {
    seq
  }
}

# Function to recursively find FASTA files with a specific pattern
find_fasta_files <- function(base_folder, pheromone_name) {
  # Create a regex pattern to match file names
  file_pattern <- paste0("^fasta_", pheromone_name, "_\\d+_nt_trimmed\\.aln$")
  
  # Search for all potential files
  all_files <- list.files(
    path = base_folder, 
    recursive = TRUE, 
    full.names = TRUE, 
    pattern = file_pattern
  )
  
  # Filter files based on parent folder name
  folder_pattern <- paste0("^", pheromone_name, "_\\d+$")
  fasta_files <- all_files[
    basename(dirname(all_files)) %in% grep(folder_pattern, basename(dirname(all_files)), value = TRUE)
  ]
  
  return(fasta_files)
}

# Function to filter out FASTA files with fewer than 2 sequences
# AND with only one unique genus_species in headers
filter_fasta_files <- function(fasta_files) {
  filtered_files <- list()
  
  for (f in fasta_files) {
    sequences <- readDNAStringSet(f)  # Read sequences
    
    # Skip file if it has fewer than 2 sequences
    if (length(sequences) < 2) {
      cat("Skipping file (fewer than 2 sequences):", f, "\n")
      next
    }
    
    # Extract genus_species from headers like "ample1_Genus_species"
    headers <- names(sequences)
    
    # Extract the Genus_species from the header using the correct regex
    species_names <- gsub("^[^_]+_([^_]+_[^_]+)$", "\\1", headers)  # Correct regex
    print(species_names)
    
    # Check if there are at least 2 unique genus_species names
    if (length(unique(species_names)) >= 2) {
      cat("File passed filter:", f, "\n")
      filtered_files <- append(filtered_files, f)  # Add file to filtered list
    } else {
      cat("Skipping file (not enough unique genus_species):", f, "\n")
    }
  }
  
  if (length(filtered_files) == 0) {
    stop("No FASTA files with at least 2 sequences and at least 2 unique genus_species found.")
  }
  
  return(filtered_files)
}

# Function to read and process all FASTA files
concatenate_alignments <- function(fasta_files, output_file) {
  alignment_list <- lapply(fasta_files, function(f) {
    readDNAStringSet(f)  # Read sequences in FASTA format
  })
  
  # Determine the length of each alignment
  gene_lengths <- sapply(alignment_list, function(x) {
    max(width(x))
  })
  
  # Collect all unique sequence headers
  all_headers <- unique(unlist(lapply(alignment_list, names)))
  
  # Initialize a list to store padded sequences
  concatenated_sequences <- setNames(vector("list", length(all_headers)), all_headers)
  
  # Loop through each alignment and pad sequences with gaps
  for (i in seq_along(alignment_list)) {
    alignment <- alignment_list[[i]]
    gene_length <- gene_lengths[i]
    
    for (header in all_headers) {
      if (header %in% names(alignment)) {
        # If sequence exists, pad to correct length
        seq <- alignment[[header]]
        padded_seq <- pad_with_gaps(seq, gene_length)
      } else {
        # If sequence is missing, fill with gaps
        padded_seq <- strrep("-", gene_length)
      }
      # Concatenate to existing sequence
      concatenated_sequences[[header]] <- paste0(concatenated_sequences[[header]], padded_seq)
    }
  }
  
  # Write concatenated sequences to output FASTA file
  final_sequences <- DNAStringSet(unlist(concatenated_sequences))
  names(final_sequences) <- names(concatenated_sequences)
  writeXStringSet(final_sequences, filepath = output_file, format = "fasta")
  
  cat("Concatenated alignment saved to:", output_file, "\n")
}

# Main Script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <base_folder> <pheromone_name>")
}

base_folder <- args[1]       # The root directory to search
pheromone_name <- args[2]    # Pheromone name like "FAD" or "FAR"

# Find all relevant FASTA files
fasta_files <- find_fasta_files(base_folder, pheromone_name)

# Filter FASTA files to retain only those with at least 2 sequences and at least 2 unique genus_species
filtered_fasta_files <- filter_fasta_files(fasta_files)

# Specify output file
output_file <- paste0("concatenated_", pheromone_name, "_alignment.fasta")

# Run the concatenation
concatenate_alignments(filtered_fasta_files, output_file)
