#!/bin/bash

#SBATCH --mem=15GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-06:00:00

# Input variables 
RES_PATH=$1
PHE=$2
COUNTER=$3

# Load IQ tree and R
module load IQ-TREE/2.3.6-gompi-2023a
module load R/4.2.1-foss-2022a

# Check number of sequences in the tree
NUM_SEQ=$(grep -c ">" $RES_PATH/fasta_${PHE}_${COUNTER}_nt_trimmed.aln)

# Conditional execution if more than 5 sequences run analyses with bootstraps and multiple threads (otherwise a single thread)
if [ "$NUM_SEQ" -gt 5 ]; then
    iqtree2 -s $RES_PATH/fasta_${PHE}_${COUNTER}_nt_trimmed.aln -m MFP -B 1000 -nt 8 --prefix $RES_PATH/${PHE}_${COUNTER}
    Rscript ./plot_tree.R --tree $RES_PATH/${PHE}_${COUNTER}.treefile --prefix $RES_PATH/${PHE}_${COUNTER}
    echo TREE FOR ${PHE} ${COUNTER} BUILT
elif [ "$NUM_SEQ" -ge 3 ]; then
    iqtree2 -s $RES_PATH/fasta_${PHE}_${COUNTER}_nt_trimmed.aln -m MFP -nt AUTO --prefix $RES_PATH/${PHE}_${COUNTER}
    Rscript ./plot_tree.R --tree $RES_PATH/${PHE}_${COUNTER}.treefile --prefix $RES_PATH/${PHE}_${COUNTER}
    echo TREE FOR ${PHE} ${COUNTER} BUILT
else
    echo "SKIPPING AS THERE ARE LESS THAN 3 SEQUENCES FOR ${PHE} ${COUNTER}."
fi

