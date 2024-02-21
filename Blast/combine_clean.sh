#!/bin/bash

#SBATCH --mem=3GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-00:10:00

RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results/$1"
cd  $RESULTS
cat blast_hits_*|grep query|head -n 1 > combine_best_hits.txt
cat blast_hits_*|grep -v query|awk ' $12 > 200' >> combine_best_hits.txt
cat blast_hits_*|grep -v query|awk ' $12 > 200 {print $1}'|sort -u > Protein_best_hits.txt
wget http://download.lepbase.org/v4/features/Heliconius_melpomene_melpomene_Hmel2.5.gff3.gz

gunzip Heliconius_melpomene_melpomene_Hmel2.5.gff3.gz
cat Protein_best_hits.txt |while read line; do grep $line Heliconius_melpomene_melpomene_Hmel2.5.gff3; done|grep mRNA|awk '{print $1"\t"$4"\t"$5}' >  location_Protein_best_hits.txt
#rm  output_chunk* blast_"$1"_hits_*
