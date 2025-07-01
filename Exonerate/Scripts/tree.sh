#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=comb


SPECIES=("${@:1:$(($#-1))}")
PHE="${@: -1}"
SCRIPT="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Scripts"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Results"

FILE="combined_${PHE}_multisp_genes.fasta"

mkdir -p $RESULTS/trees/${PHE}
cd $RESULTS

if [ -f "$RESULTS/trees/${PHE}/$FILE" ]; then
    rm $RESULTS/trees/${PHE}/"$FILE"
fi
touch $RESULTS/trees/${PHE}/"$FILE"

# Combine all alignment
for i in "${SPECIES[@]}"
   do cat $i/${PHE}/unique_${i}_${PHE}_combined_genes.fasta >> $RESULTS/trees/${PHE}/$FILE
done

# Align with muscle
module load MUSCLE/3.8.1551-GCC-9.3.0
pwd
muscle -in  $RESULTS/trees/${PHE}/$FILE -out  $RESULTS/trees/${PHE}/muscle_aln
module load MAFFT/7.505-GCC-11.3.0-with-extensions
mafft --maxiterate 1000 --localpair  $RESULTS/trees/${PHE}/muscle_aln > $RESULTS/trees/${PHE}/aligned_"$FILE"
module load Biopython/1.83-foss-2023a
cat $RESULTS/trees/${PHE}/aligned_"$FILE"|python3 $SCRIPT/filter_too_short_long_sequence_excluding_gaps.py -o $RESULTS/trees/${PHE}/aligned_"$FILE"_filtered -s kept_sequences -r 0.3 -m 1.5
module purge
module load R/4.2.1-foss-2022a
export PATH=~/local/bin:$PATH
Rscript $SCRIPT/msa_plot.R  $RESULTS/trees/${PHE}/aligned_"$FILE"_filtered $RESULTS/trees/${PHE}/alignment_msa.html

# trim the alignment
module purge
module load trimAl/1.4.1-GCC-9.3.0
trimal -in $RESULTS/trees/${PHE}/aligned_"$FILE"_filtered  -gt 0.7 -cons 20 -out $RESULTS/trees/${PHE}/trimal_aligned_"$FILE"

module purge
module load R/4.2.1-foss-2022a
export PATH=~/local/bin:$PATH
Rscript $SCRIPT/msa_plot.R  $RESULTS/trees/${PHE}/trimal_aligned_"$FILE" $RESULTS/trees/${PHE}/trimal_alignment_msa.html

# Infer tree
module load IQ-TREE/2.3.6-gompi-2023a
iqtree2 -s $RESULTS/trees/${PHE}/trimal_aligned_"$FILE" -m MFP -B 1000 --prefix $RESULTS/trees/${PHE}/tree_${PHE} -redo
	
# Plot tree with ggtree
module load R/4.2.1-foss-2022a
cp $SCRIPT/species_colors.csv ./
Rscript $SCRIPT/plot_tree.R --tree $RESULTS/trees/${PHE}/tree_${PHE}.treefile --prefix $RESULTS/trees/${PHE}/tree_${PHE}



