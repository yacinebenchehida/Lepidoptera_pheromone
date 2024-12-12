#!/usr/bin/env Rscript

# Load required libraries
library(ggtree)
library(ape)
library(argparse)
library(dplyr)

# Create argument parser
parser <- ArgumentParser(description = "Plot unrooted tree from a Newick file")
parser$add_argument("-t", "--tree", type = "character", required = TRUE, help = "Path to the tree file (Newick format)")
parser$add_argument("-p", "--prefix", type = "character", required = TRUE, help = "Output path and prefix for the tree plots")
parser$add_argument("--color_file", type = "character", default = "species_colors.csv", help = "Path to save/load the consistent species-color mapping file")

# Parse arguments
args <- parser$parse_args()

# Load tree from Newick file
tree <- read.tree(args$tree)

# Extract species labels (based on the last part of the label after the last underscore)
species_labels <- gsub("^[^_]+_([^_]+_[^_]+)$", "\\1", tree$tip.label)  # Extract genus_species from the label
tree$tip.label <- gsub("_[^_]+_[^_]+$", "", tree$tip.label)  # Remove the genus_species part from the label for the tree plot

# Convert tree to ggtree object for easier handling
tree2 <- ggtree(tree, layout = "unrooted")

# Create a factor for species to ensure consistent coloring
species_colors <- factor(species_labels)

# Check if the color mapping file exists
if (file.exists(args$color_file)) {
  # Load the existing species-color mapping
  species_color_map <- read.csv(args$color_file)
  
  # Ensure the CSV contains both species and their corresponding colors
  if(!all(c("Species", "Color") %in% colnames(species_color_map))) {
    stop("The CSV file must contain 'Species' and 'Color' columns.")
  }

  # Create a named vector mapping species to colors
  color_palette <- setNames(species_color_map$Color, species_color_map$Species)
  
} else {
  # If the color mapping file does not exist, create a new color palette
  color_palette <- rainbow(length(unique(species_colors)))
  
  # Create a data frame with species and their associated color
  species_color_map <- data.frame(Species = unique(species_colors),
                                  Color = color_palette, stringsAsFactors = FALSE)
  
  # Save the species-color mapping to a file for future use
  write.csv(species_color_map, args$color_file, row.names = FALSE)
}

# Add species labels to the tree$data dataframe (only for tips)
tree2$data$species <- NA
tree2$data$species[tree2$data$isTip] <- species_labels  # Only assign species labels to tips

# Merge species labels with the color mapping (to ensure we have the right colors for each species)
tree2$data <- left_join(tree2$data, species_color_map, by = c("species" = "Species"))

# Plot tree using ggtree and map species to tip labels (with unrooted layout)
p <- ggtree(tree2$data, layout = "unrooted") +  
  geom_tiplab(aes(label = label, color = species), size = 2) +  # Align labels and color by species
  scale_color_manual(values = color_palette)  

# Define output file paths
pdf_file <- paste0(args$prefix, ".pdf")
png_file <- paste0(args$prefix, ".png")

# Save the tree plots
ggsave(pdf_file, plot = p, device = "pdf", width = 16, height = 12)
ggsave(png_file, plot = p, device = "png", width = 16, height = 12)

cat("Tree plot saved to", pdf_file, "\n")
