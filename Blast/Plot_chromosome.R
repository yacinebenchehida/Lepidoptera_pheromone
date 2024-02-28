#!/bin/bash

#SBATCH --mem=3GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-00:10:00


module load R/4.2.1-foss-2022a

Rscript --vanilla Plot_chromosome.R $1
[ybc502@login2[viking2] Script]$ cat Plot_chromosome.R
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
Karytype = read.table(paste("../Inputs/",SPECIES,"/scaffold_size_information.txt",sep=""),sep = "\t")
#Karytype$V1 = gsub(x = Karytype$V1, pattern = "Hmel2(\\d{2})(.+)",replacement = "\\1",perl = TRUE)
colnames(Karytype) = c("Chr","Start","End")
print("genome size file read")

#########################
# Read annotation files #
#########################
anno_files <- lapply(system(paste("readlink -f ../Results/",SPECIES,"/*/location*",sep=""),intern = T), read.table, header = FALSE, sep = "\t")
anno_files = as.data.frame(do.call(rbind,anno_files))
#anno_files$V2 =  gsub(x = anno_files$V2, pattern = "Hmel2(\\d{2})(.+)",replacement = "\\1",perl = TRUE)
colnames(anno_files) = c("Gene","Chr","Start","End","Type","BitScore")
anno_files$BitScore = as.numeric(as.character(anno_files$BitScore))
print("annotation files read")

################################
# Adjust Bitscore for plotting #
################################
anno_files[anno_files$Type=="FAD",6] = anno_files[anno_files$Type=="FAD",6] - 150
anno_files[anno_files$Type=="FAR",6] = anno_files[anno_files$Type=="FAR",6] - 350

######################################################
# Add genome size information to the annotation file #
######################################################
anno_files$Size = 0

for (i in Karytype$Chr){
  anno_files[which(grepl(i, anno_files$Chr)),7] = Karytype[Karytype$Chr==i,3]
}
print("information about genome size added to annotation files")

##############################################
# Remove genes falling in very small contigs #
##############################################
anno_files = anno_files[anno_files$Size > 1000000,]
Karytype = Karytype[Karytype$End > 1000000,]
print(Karytype)

##############################################
# Rename chromosome as 1,2,3,4,5,6,7,8,9.... #
##############################################
Karytype$NEW_name = 0
for(i in 1:dim(Karytype)[1]){
  Karytype[i,4] = i
}

Karytype$NEW_name <- factor(Karytype$NEW_name, levels = sort(unique(as.numeric(Karytype$NEW_name))))

for (i in Karytype$Chr){
  anno_files[which(grepl(i, anno_files$Chr)),2] = Karytype[Karytype$Chr==i,4]
}

anno_files$Chr <- factor(anno_files$Chr, levels = sort(unique(as.numeric(anno_files$Chr))))
Karytype = Karytype[,-c(1)]
colnames(Karytype) = c("Start","End","Chr")
print(anno_files)

########
# Plot #
########
print("start plotting")
pdf(paste("../Results/",SPECIES,"/",SPECIES,"_FAD_FAR_Location.pdf",sep=""),7,dim(Karytype)[1]/4)
ggplot(data = Karytype, aes(xmin = (0-200000)/1000000, xmax = (End+50000)/1000000,ymin = 0, ymax = 0.5)) +
  ggchicklet:::geom_rrect(fill = "white", colour = "black",) + guides(color = 'none') + theme_bw() + 
  facet_grid(Chr~., switch= "y", drop =TRUE,labeller = label_parsed, scales = "free", space = "free") +
  theme(axis.ticks.y = element_blank(), strip.placement = "outside", axis.text.y = element_blank(), panel.spacing = unit(0, "lines")) +
  xlab("Position (Mb)") + ylab("Chromosome") + 
  geom_rect(data = anno_files, aes(xmin = (Start-25000)/1000000, xmax = (End+25000)/1000000, ymin = 0, ymax = 0.5, fill = Type, alpha = BitScore)) +
  ggtitle(SPECIES) + theme(plot.title = element_text(hjust = 0.5))
dev.off()
