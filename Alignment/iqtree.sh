#!/bin/bash

#SBATCH --mem=15GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-06:00:00

RES_PATH=$1
PHE=$2
COUNTER=$3

module load IQ-TREE/2.3.6-gompi-2023a

iqtree2 -s $RES_PATH/fasta_${PHE}_${COUNTER}_nt_trimmed.aln -m MFP -B 1000 -nt 8 --prefix $RES_PATH/${PHE}_${COUNTER}
echo TREE FOR ${PHE} ${COUNTER} BUILT
