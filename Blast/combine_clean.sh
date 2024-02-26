#!/bin/bash

#SBATCH --mem=3GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-00:10:00

RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results/$1"
PHE=$1
BITSCORE=$2
cd  $RESULTS
cat blast_hits_*|grep query|head -n 1 > combine_best_hits.txt
cat blast_hits_*|grep -v query|awk -v biscore=$BITSCORE '$12 > biscore' >> combine_best_hits.txt
cat blast_hits_*|grep -v query|awk -v biscore=$BITSCORE '$12 > biscore {print $1}'|sort -u > Protein_best_hits.txt

# Annotation file downloaded
wget http://download.lepbase.org/v4/features/Heliconius_melpomene_melpomene_Hmel2.5.gff3.gz 
echo ANNOTATION FILE DOWNLOADED

gunzip Heliconius_melpomene_melpomene_Hmel2.5.gff3.gz
echo ANNOTATION FILE DECOMPRESSED

#cat Protein_best_hits.txt |while read line; do grep $line Heliconius_melpomene_melpomene_Hmel2.5.gff3; done|grep mRNA|awk '{print $1"\t"$4"\t"$5}' >  location_Protein_best_hits.txt
cat Protein_best_hits.txt | while read line; do paste <(echo "$line") <(grep "$line" Heliconius_melpomene_melpomene_Hmel2.5.gff3 | grep mRNA) | awk -v PHE="$PHE" '{print $1"\t"$2"\t"$5"\t"$6"\t"PHE"\t"$12 }'; done > location_Protein_best_hits.txt
echo LOCATION OF $1 written
#rm  output_chunk* blast_"$1"_hits_*
rm output_chunk*

python ../../Script/extract_fasta.py Protein_best_hits.txt Heliconius_melpomene_proteins.fa fasta_"$1"_protein.fa
