#################
# Load packages #
#################
library(ggchicklet)
library(ggplot2)

#########################
# Read annotation files #
#########################
dat = commandArgs(trailingOnly=TRUE)
SPECIES=dat
anno_files <- read.table("gene_position_to_plot.txt",sep="\t")
colnames(anno_files) = c("Chr","Start","End")
print("annotation files read")

##############################################
# Function to keep non overlapping intervals #
##############################################
merge_intervals <- function(df) {
  colnames(df) <- c("chr", "start", "end")
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  
  result <- data.frame()
  
  for (chrom in unique(df$chr)) {
    sub <- df[df$chr == chrom, ]
    sub <- sub[order(sub$start), ]
    
    if (nrow(sub) == 1) {
      # Only one interval, no merging needed
      result <- rbind(result, sub)
      next
    }
    
    merged <- list()
    current_start <- sub$start[1]
    current_end <- sub$end[1]
    
    for (i in 2:nrow(sub)) {
      s <- sub$start[i]
      e <- sub$end[i]
      
      if (s <= current_end) {
        current_end <- max(current_end, e)
      } else {
        merged[[length(merged) + 1]] <- c(chrom, current_start, current_end)
        current_start <- s
        current_end <- e
      }
    }
    merged[[length(merged) + 1]] <- c(chrom, current_start, current_end)
    
    merged_df <- do.call(rbind, merged)
    colnames(merged_df) <- c("chr", "start", "end")
    result <- rbind(result, merged_df)
  }
  
  result <- as.data.frame(result)
  result$start <- as.numeric(result$start)
  result$end <- as.numeric(result$end)
  
  return(result)
}
  
anno_files <- merge_intervals(anno_files)
colnames(anno_files) <- c("Chr","Start","End")
write.table(x = anno_files, "annotation_merged.txt", quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)
print(anno_files)