# Load required libraries
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

# Function to read and process all FASTA files
concatenate_alignments <- function(fasta_files, output_file) {
  # Step 1: Read all files into a list
  alignment_list <- lapply(fasta_files, function(f) {
    readDNAStringSet(f)  # Reads sequences in FASTA format
  })
  
  # Step 2: Determine the length of each alignment
  gene_lengths <- sapply(alignment_list, function(x) {
    max(width(x))
  })
  total_length <- sum(gene_lengths)
  
  # Step 3: Collect all unique sequence headers
  all_headers <- unique(unlist(lapply(alignment_list, names)))
  
  # Step 4: Initialize a list to store padded sequences
  concatenated_sequences <- setNames(vector("list", length(all_headers)), all_headers)
  
  # Step 5: Loop through each alignment and pad sequences with gaps
  for (i in seq_along(alignment_list)) {
    alignment <- alignment_list[[i]]
    gene_length <- gene_lengths[i]
    
    for (header in all_headers) {
      if (header %in% names(alignment)) {
        # If the sequence exists for this gene, pad to correct length
        seq <- alignment[[header]]
        padded_seq <- pad_with_gaps(seq, gene_length)
      } else {
        # If sequence is missing, fill with gaps
        padded_seq <- strrep("-", gene_length)
      }
      # Concatenate to existing sequence
      if (is.null(concatenated_sequences[[header]])) {
        concatenated_sequences[[header]] <- padded_seq
      } else {
        concatenated_sequences[[header]] <- paste0(concatenated_sequences[[header]], padded_seq)
      }
    }
  }
  
  # Step 6: Write concatenated sequences to output FASTA file
  final_sequences <- DNAStringSet(unlist(concatenated_sequences))
  names(final_sequences) <- names(concatenated_sequences)
  writeXStringSet(final_sequences, filepath = output_file, format = "fasta")
  
  cat("Concatenated alignment saved to:", output_file, "\n")
}

# Main Script
# Set base folder and pheromone name as arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript script.R <base_folder> <pheromone_name>")
}

base_folder <- args[1]       # The root directory to search
pheromone_name <- args[2]    # Pheromone name like "FAD" or "FAR"

# Find all relevant FASTA files
fasta_files <- find_fasta_files(base_folder, pheromone_name)

# Check if files were found
if (length(fasta_files) == 0) {
  stop("No matching FASTA files found.")
}

# Specify output file
output_file <- paste0("concatenated_", pheromone_name, "_alignment.fasta")

# Run the concatenation
concatenate_alignments(fasta_files, output_file)
