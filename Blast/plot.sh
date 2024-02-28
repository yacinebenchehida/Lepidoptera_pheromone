#!/bin/bash

#SBATCH --mem=3GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-00:10:00


module load R/4.2.1-foss-2022a

Rscript --vanilla Plot_chromosome.R $1
