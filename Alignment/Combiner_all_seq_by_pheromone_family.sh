#!/bin/bash

#SBATCH --mem=15GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-06:00:00

# Input variables 
CHEMIN=$1
PHE=$2

# Load IQ tree and R
module load IQ-TREE/2.3.6-gompi-2023a
module load R/4.2.1-foss-2022a

mkdir -p ${CHEMIN}/All_combined
RES_PATH=${CHEMIN}/All_combined
echo $RES_PATH

# Combine all data per pheromone
Rscript ./combine_all_seq_by_pheromone_family.R ${CHEMIN}/${PHE} ${PHE}
mv concatenated_${PHE}_alignment.fasta $RES_PATH
echo FASTA COMBINED

# Get Phylogenetic tree
iqtree2 -s $RES_PATH/concatenated_${PHE}_alignment.fasta -m MFP -B 1000 -nt 8 --prefix $RES_PATH/All_${PHE}
echo IQTREE DONE

# Plot Tree
Rscript ./plot_tree.R --tree $RES_PATH/All_${PHE}.treefile --prefix $RES_PATH/All_${PHE}
echo FINAL COMBINED TREE GENERATED
