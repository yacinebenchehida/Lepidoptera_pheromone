# Load packages
library(ggchicklet)
library(ggplot2)

# Read genome size files
Heliconius_Karytype = read.table("../Inputs/scaffold_size_information.txt",sep = "\t")
Heliconius_Karytype$V1 = gsub(x = Heliconius_Karytype$V1, pattern = "Hmel2(\\d{2})(.+)",replacement = "\\1",perl = TRUE)
colnames(Heliconius_Karytype) = c("Chr","Start","End")
print("genome size file read")

# Read annotation files
anno_files <- lapply(system("readlink -f ../Results/*/location*",intern = T), read.table, header = FALSE, sep = "\t")
anno_files = as.data.frame(do.call(rbind,anno_files))
anno_files$V2 =  gsub(x = anno_files$V2, pattern = "Hmel2(\\d{2})(.+)",replacement = "\\1",perl = TRUE)
colnames(anno_files) = c("Gene","Chr","Start","End","Type")
anno_files$Chr = as.factor(anno_files$Chr)
print("annotation files read")

# Add genome size information to the annotation file
anno_files$Size = 0

for (i in Heliconius_Karytype$Chr){
  anno_files[which(grepl(i, anno_files$Chr)),6] = Heliconius_Karytype[Heliconius_Karytype$Chr==i,3]
}
print("information about genome size added to annotation files")


# Remove genes falling in very small contigs, 
anno_files = anno_files[anno_files$Size > 0,]

# Plot
print("start plotting")
pdf("FAD_FAR_Location.pdf",7,7)
ggplot(data = Heliconius_Karytype, aes(xmin = (0-200000)/1000000, xmax = (End)/1000000,ymin = 0, ymax = 0.5)) +
  ggchicklet:::geom_rrect(fill = "white", colour = "black",) + guides(color = 'none') + theme_bw() + 
  facet_grid(Chr~., switch= "y", drop =TRUE,labeller = label_parsed, scales = "free", space = "free") +
  theme(axis.ticks.y = element_blank(), strip.placement = "outside", axis.text.y = element_blank(), panel.spacing = unit(0, "lines")) +
  xlab("Position (Mb)") + ylab("Chromosome") +
  geom_rect(data = anno_files, aes(xmin = (Start-30000)/1000000, xmax = (End+30000)/1000000, ymin = 0, ymax = 0.5, fill = Type))
dev.off()
