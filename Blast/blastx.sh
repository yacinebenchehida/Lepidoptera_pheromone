#!/bin/bash

#SBATCH --mem=5GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-2:00:00
#SBATCH --job-name=blastx

RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results/$1"
python3 ./blastx_multi_fasta.py $RESULTS/output_chunk_${SLURM_ARRAY_TASK_ID}.fasta ${SLURM_ARRAY_TASK_ID} $RESULTS/"$1"_db
