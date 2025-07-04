#################
# Load packages #
#################
library(ggchicklet)
library(ggplot2)

##########################
# Read genome size files #
##########################
dat = commandArgs(trailingOnly=TRUE)
SPECIES=dat
Karytype = read.table("scaffold_size_information.txt",sep = "\t")
#Karytype$V1 = gsub(x = Karytype$V1, pattern = "Hmel2(\\d{2})(.+)",replacement = "\\1",perl = TRUE)
colnames(Karytype) = c("Chr","Start","End")
print("genome size file read")

#########################
# Read annotation files #
#########################
anno_files <- read.table("genes_2_plot.txt",sep="\t")
#anno_files$V2 =  gsub(x = anno_files$V2, pattern = "Hmel2(\\d{2})(.+)",replacement = "\\1",perl = TRUE)
colnames(anno_files) = c("Chr","Start","End")
print("annotation files read")

######################################################
# Add genome size information to the annotation file #
######################################################
anno_files$Size = 0

for (i in Karytype$Chr){
  anno_files[which(grepl(paste(i,"$",sep=""), anno_files$Chr)),4] = Karytype[Karytype$Chr==i,3]
}
print("information about genome size added to annotation files")

##############################################
# Remove genes falling in very small contigs #
##############################################
anno_files = anno_files[anno_files$Size > 1000000,]
Karytype = Karytype[Karytype$End > 1000000,]

##############################################
# Rename chromosome as 1,2,3,4,5,6,7,8,9.... #
##############################################
# Try to coerce chromosome names to numeric (NA if not numeric)
numeric_chr <- suppressWarnings(as.numeric(as.character(Karytype$Chr)))
is_numeric <- !is.na(numeric_chr)

# Numeric chromosome IDs kept as is
numeric_ids <- numeric_chr[is_numeric]

# Non-numeric chromosomes get assigned numbers starting after max numeric chromosome
non_numeric_chr <- Karytype$Chr[!is_numeric]
start_num <- ifelse(length(numeric_ids) > 0, max(numeric_ids), 0)
non_numeric_ids <- seq(start_num + 1, length.out = length(non_numeric_chr))

# Prepare NEW_name vector of characters
NEW_name <- character(nrow(Karytype))
NEW_name[is_numeric] <- as.character(numeric_ids)
NEW_name[!is_numeric] <- as.character(non_numeric_ids)

Karytype$NEW_name <- factor(NEW_name, levels = sort(as.numeric(unique(NEW_name))))

# Update annotation file chromosome names accordingly
for (i in seq_len(nrow(Karytype))) {
  old_chr <- Karytype$Chr[i]
  new_chr <- as.character(Karytype$NEW_name[i])
  anno_files$Chr[anno_files$Chr == old_chr] <- new_chr
}

anno_files$Chr <- factor(anno_files$Chr, levels = levels(Karytype$NEW_name))

Karytype = Karytype[,-c(1)]
colnames(Karytype) = c("Start","End","Chr")

########
# Plot #
########
print("start plotting")
pdf(paste(SPECIES,"_Location.pdf",sep=""),7,dim(Karytype)[1]/4)
ggplot(data = Karytype, aes(xmin = (0-200000)/1000000, xmax = (End+50000)/1000000,ymin = 0, ymax = 0.5)) +
  ggchicklet:::geom_rrect(fill = "white", colour = "black",) + guides(color = 'none') + theme_bw() + 
  facet_grid(Chr~., switch= "y", drop =TRUE,labeller = label_parsed, scales = "free", space = "free") +
  theme(axis.ticks.y = element_blank(), strip.placement = "outside", axis.text.y = element_blank(), panel.spacing = unit(0, "lines")) +
  xlab("Position (Mb)") + ylab("Chromosome") + 
  geom_rect(data = anno_files, aes(xmin = (Start/1000000)-5000/1000000, xmax = (End/1000000)+5000/1000000, ymin = 0, ymax = 0.5, fill = "black")) +
  ggtitle(SPECIES) + theme(plot.title = element_text(hjust = 0.5))
dev.off()
