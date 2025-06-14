#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-0:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=comb

SCRIPT="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Scripts"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Results"
SPECIES=$@

FILE="combined_FAD_multisp_genes.fasta"

cd $RESULTS
if [ -f "$FILE" ]; then
    rm "$FILE"
fi
touch "$FILE"

# Combine all alignment

for i in $SPECIES
    do cat $i/FAD/trimal_species_name_${i}_FAD_combined_genes.fasta >> $FILE
done

# Align with muscle
module load MUSCLE/3.8.1551-GCC-9.3.0
muscle -in $FILE -out aligned_"$FILE"
module load Biopython/1.83-foss-2023a
cat aligned_"$FILE"|python3 $SCRIPT/filter_too_short_long_sequence.py -o aligned_"$FILE"_filtered -s kept_sequences -r 0.4 -m 2.0
module purge
module load R/4.2.1-foss-2022a
export PATH=~/local/bin:$PATH
Rscript $SCRIPT/msa_plot.R  aligned_"$FILE"_filtered alignment_msa.html

# trim the alignment
module purge
module load trimAl/1.4.1-GCC-9.3.0
trimal -in aligned_"$FILE"_filtered  -gt 0.7 -cons 20 -out trimal_aligned_"$FILE"

module purge
module load R/4.2.1-foss-2022a
export PATH=~/local/bin:$PATH
Rscript $SCRIPT/msa_plot.R  trimal_aligned_"$FILE" trimal_alignment_msa.html

# Infer tree
module load IQ-TREE/2.3.6-gompi-2023a
iqtree2 -s trimal_aligned_"$FILE" -m MFP -B 1000 --prefix tree_FAD -redo
	
# Plot tree with ggtree
module load R/4.2.1-foss-2022a
cp $SCRIPT/species_colors.csv ./
Rscript $SCRIPT/plot_tree.R --tree tree_FAD.treefile --prefix tree_FAD



