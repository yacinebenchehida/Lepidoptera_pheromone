#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=comb


SCRIPT="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Scripts"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Results/$1"
DB="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results/Heliconius_melpomene/FAD/FAD_db"
DATA="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Inputs/$1"
FASTA_REF="${1}_genome.fa"
REF_GENOME="$DATA/$FASTA_REF"
echo $REF_GENOME

cd $RESULTS/$2
ls

# Merge results of each gene inferred by exonerate for each species
FILE1="${1}_${2}_combined_genes.fasta"

if [ -f "$FILE1" ]; then
    rm "$FILE1"
fi
touch "$FILE1"

for i in output_chunk_*; do
	cat $i/*_cds >> "$FILE1"
done


# add species name
cat "$FILE1" | awk -v species="_${1}" '/^>/ {split($1, a, " "); print a[1] species; next} {print}' > unique_sp_name_cds_"$FILE1"
module load Biopython/1.83-foss-2023a
python3 $SCRIPT/unique_sequence_1_pick.py unique_sp_name_cds_"$FILE1" > unique_species_name_cds_"$FILE1"
cat unique_species_name_cds_"$FILE1"|python3 $SCRIPT/filter_too_short_long_sequence.py -o "$FILE1"_filtered -s kept_sequences -r 0.90 -m 2.0
python3 $SCRIPT/extract_longuest_fasta_based_on_overlapping_interval.py "$FILE1"_filtered > unique_interval_"$FILE1"_filtered

# Check cds sequenes by reblasting them against the reference database
module load BLAST+/2.14.1-gompi-2023a
blastp -query unique_interval_"$FILE1"_filtered -db $DB -outfmt 6 -max_target_seqs 1  -evalue 1e-20 > blast_results
#python3 $SCRIPT/filter_blast_by_overlap_bitscore.py blast_results > non_overlaping_blast_results.txt
echo BLASTP DONE

# Plot genes on ideogram
python3 $SCRIPT/scaffold_size.py $REF_GENOME|awk '$3 > 1000000' > scaffold_size_information.txt

module purge
module load R/4.2.1-foss-2022a
#cat non_overlaping_blast_results.txt| cut -f1-3 -d _ |perl -pe 's/_/\t/g' > genes_2_plot.txt
 
#cp $SCRIPT/Plot_chromosome.R ./
#Rscript ./Plot_chromosome.R ${1}
#rm Plot_chromosome.R
#echo GOOD GENE PLOTTED

# Extract genes that blast properly to the FAD or FAD database
cat blast_results|awk '{print $1}' > matches.txt
awk 'NR==FNR {ids[$1]; next} /^>/ {header=$0; id=substr($1,2); keep=(id in ids)} keep' matches.txt  unique_species_name_cds_"$FILE1" > checked_unique_species_name_cds_"$FILE1"
module load Biopython/1.83-foss-2023a
python3 $SCRIPT/unique_sequence_1_pick.py checked_unique_species_name_cds_"$FILE1" > checked_unique_sp_name_cds_"$FILE1"

# align sequences with muscle
module load MUSCLE/3.8.1551-GCC-9.3.0
muscle -in unique_interval_"$FILE1"_filtered -out aligned_unique_species_name_cds_"$FILE1"
echo MUSCLE ALIGNMENT DONE

# plot alignment
module purge
module load R/4.2.1-foss-2022a
export PATH=~/local/bin:$PATH
Rscript $SCRIPT/msa_plot.R aligned_unique_species_name_cds_"$FILE1" alignment_msa_${1}_${2}.html
echo MSA PLOT BEFORE TRIMAL GENERATED

# Trim sequence with trimal
module purge
module load trimAl/1.4.1-GCC-9.3.0

trimal -in aligned_unique_species_name_cds_"$FILE1" -gt 0.5 -cons 20 -out trimal  # -gt 0.7 more than 1-0.7=0.3 so removes sites with more than 30% of gaps. -cons 30  
#trimal -in aligned_unique_species_name_cds_"$FILE1" -gappyout -out trimal
cat  trimal | awk -v species="_${1}" '/^>/ {split($1, a, " "); print a[1] species; next} {print}' > trimal_species_name_"$FILE1"
echo TRIMAL DONE

module purge
module load Biopython/1.83-foss-2023a 
cat trimal_species_name_"$FILE1"|python3 $SCRIPT/filter_too_short_long_sequence_excluding_gaps.py -o trimal_"$FILE1"_filtered -s kept_sequences -r 0.9 -m 2.0


module purge
module load R/4.2.1-foss-2022a
export PATH=~/local/bin:$PATH
Rscript $SCRIPT/msa_plot.R trimal_"$FILE1"_filtered trimal_alignment_msa_${1}_${2}.html
echo MSA PLOT AFTER TRIMAL GENERATED