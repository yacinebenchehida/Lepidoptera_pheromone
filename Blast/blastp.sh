#!/bin/bash

#SBATCH --mem=5GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-2:00:00

python3 ./blastp_multi_fasta.py ../Results/output_chunk_${SLURM_ARRAY_TASK_ID}.fasta ${SLURM_ARRAY_TASK_ID}
