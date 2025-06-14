#!/usr/bin/env Rscript

# Load required libraries
library(ggtree)
library(ape)
library(argparse)
library(dplyr)
library(ggplot2)
library(ggnewscale)  # Import ggnewscale for handling multiple color scales

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
species_labels <- gsub(".*_([^_]+_[^_]+)$", "\\1", tree$tip.label)  # Extract genus_species from the label
tree$tip.label <- gsub("_[^_]+_[^_]+$", "", tree$tip.label)  # Simplify tip labels
print(tree$tip.label)

plotting <- function(methode){
# Convert tree to ggtree object
tree2 <- ggtree(tree, layout = methode)

if (methode %in% c("circular", "daylight","rectangular")) {
  tree2 <- ggtree(tree, layout = methode, branch.length = 'none')
}

# Assign species labels to tips
tree2$data$species <- NA
tree2$data$species[tree2$data$isTip] <- species_labels

# Check if bootstrap values exist
bootstrap_exists <- !is.null(tree$node.label) && any(tree$node.label != "")

# Assign bootstrap values to internal nodes
tree2$data$bootstrap <- NA
if (bootstrap_exists) {
  tree2$data$bootstrap[!tree2$data$isTip] <- as.numeric(tree$node.label)
}

# Handle species color mapping
species_colors <- factor(species_labels)
if (file.exists(args$color_file)) {
  species_color_map <- read.csv(args$color_file)
  if (!all(c("Species", "Color") %in% colnames(species_color_map))) {
    stop("The CSV file must contain 'Species' and 'Color' columns.")
  }
  color_palette <- setNames(species_color_map$Color, species_color_map$Species)
} else {
  color_palette <- rainbow(length(unique(species_colors)))
  species_color_map <- data.frame(Species = unique(species_colors),
                                  Color = color_palette, stringsAsFactors = FALSE)
  write.csv(species_color_map, args$color_file, row.names = FALSE)
}

# Merge species labels with colors
tree2$data <- left_join(tree2$data, species_color_map, by = c("species" = "Species"))

# Start the plot
p <- ggtree(tree2$data, layout = methode) +  
  geom_tiplab(aes(label = label, color = species), size = 2) +
  scale_color_manual(values = color_palette) +
  theme(
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.4, "cm"),
    legend.spacing.y = unit(0.1, "cm"),
    plot.margin = margin(1, 1, 1, 1, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank()
  )


# Add a new color scale for bootstrap values (continuous scale)
if (bootstrap_exists) {
  p <- p +
    ggnewscale::new_scale_fill() +  # Reset color scale
    geom_point2(aes(subset = (!isTip & !is.na(bootstrap)), x = x, y = y, fill = bootstrap),
                shape = 21, size = 2, color = "black") +  # Circles for bootstrap
    scale_fill_gradientn(colours = c("green","yellow","red"),limits = c(0, 100))  #+  # Bootstrap color gradient
     #geom_nodelab(aes(label = bootstrap), size = 3, color = "black",nudge_y = 0, nudge_x = 0.012) 
  }

  # Save the output
pdf_file <- paste0(args$prefix, "_", methode, ".pdf")
png_file <- paste0(args$prefix, "_", methode, ".png")
ggsave(pdf_file, plot = p, device = "pdf", width = 10, height = 8)
ggsave(png_file, plot = p, device = "png", width = 10, height = 8)

}


for (i in c("circular","equal_angle","slanted","daylight","rectangular")){
  plotting(i)  
}
