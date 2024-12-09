#!/bin/bash

#SBATCH --mem=5GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-00:20:00

module load Biopython/1.81-foss-2022b
python3 ./calculate_identity.py $1
echo IDENTITY CALCULATED
module purge
module load Cython/3.0.8-GCCcore-12.2.0
python -c "import filter_orthologs; filter_orthologs.main()" /mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Alignment/Results/FAD/Alignments  300 0.20 100 > Candidate_orthogs.txt
