#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-2:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=master

###########################
# 0 - Load useful modules #
###########################
module load Exonerate/2.4.0-GCC-11.3.0
module load SeqKit/2.3.1
module load Biopython/1.81-foss-2022b
module load BEDTools/2.31.0-GCC-12.3.0


#################################################
# 1 - Set reference genomes and important paths #
#################################################
DATA="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Inputs/$1"
FASTA_REF="${1}_genome.fa"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Results/$1"
Proteins="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Input/${2}/unique_db_${2}.fasta"
SCRIPT="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Scripts"
REF_GENOME="$DATA/$FASTA_REF"
echo $REF_GENOME

######################################################################
# 2 - divide FAD/FAR data base into chuncks to speed up the pipeline #
######################################################################
#python3 ./divide_fasta.py $Proteins

#####################
# 3 - Run EXONERATE #
#####################
SP=$(echo $1|cut -c 3)
jobname="${SP}exo"
array_number=$(ls /mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Input/${2}/output_chunk_*|wc -l)
sbatch --job-name="$jobname" --array=1-$array_number ./exonerate.sh $1 $2

########################################
# 4 - Combine exonerate gene sequences #
########################################
running_jobs_exonerate=$(squeue|grep $(whoami)| grep -P $jobname | awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
echo $running_jobs_exonerate
sbatch --dependency=aftercorr:$running_jobs_exonerate ./combine_results.sh $1 $2

