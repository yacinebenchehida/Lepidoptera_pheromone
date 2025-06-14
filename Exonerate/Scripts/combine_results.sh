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
	cat $i/*_tmp_genes >> "$FILE1"
done

# Get overlapping genes inferred by exonerate
grep ">" "$FILE1" | perl -pe 's/>//g'| awk -F "_" '{print $1"\t"$3"\t"$4}' > gene_position_to_plot.txt
module load Biopython/1.83-foss-2023a
python3 $SCRIPT/scaffold_size.py $REF_GENOME|awk '$3 > 1000000' > scaffold_size_information.txt
cp $SCRIPT/interval_parse.R ./
module purge
module load R/4.2.1-foss-2022a
Rscript ./interval_parse.R ${1}
rm interval_parse.R
# => return a file called annotation_merged.txt that combines overlapping genes. 

# Get a fasta file only with one fasta for each non overlapping region
module purge
module load BEDTools/2.31.0-GCC-12.3.0
bedtools getfasta -fi $REF_GENOME -bed annotation_merged.txt -fo unique_"$FILE1"
echo FASTA FILE WITH NON OVERLAPPING GENES GENERATED

# Get coding region using augustus
module purge 
module load AUGUSTUS/3.5.0-foss-2022a
augustus --progress=true --strand=both --species=heliconius_melpomene1 unique_"$FILE1" > cds_"$FILE1".gff
#augustus --progress=true --strand=both --species=fly unique_"$FILE1" > cds_"$FILE1".gff
getAnnoFasta.pl cds_"$FILE1".gff
echo AUGUSTUS DONE

# add species name
cat cds_"$FILE1".aa | awk -v species="_${1}" '/^>/ {split($1, a, " "); print a[1] species; next} {print}' > unique_species_name_cds_"$FILE1"
echo SPECIES NAME ADDED

# Check cds sequenes by reblasting them against the reference database
module load BLAST+/2.14.1-gompi-2023a
blastp -query unique_species_name_cds_"$FILE1" -db $DB -outfmt 6 -max_target_seqs 1  -evalue 1e-30 > blast_results
echo BLASTP DONE

# Plot final genes
module purge
module load R/4.2.1-foss-2022a
cat blast_results | awk '$1 > 100 {print $1}' > matches.txt
cat matches.txt | cut -f1 -d_ | cut -f1 -d"." | while read line; do     grep -P "$line$" cds_"$FILE1".gff; done | grep -P "\tgene\t" | awk 'BEGIN{OFS="\t"}
{
    split($1, a, /[:\-]/);
    scaffold = a[1];
    orig_start = a[2];
    orig_end = a[3];
    gene_start = $4;
    gene_len = $5;
    new_start = orig_start + gene_start - 1;
    new_end = orig_start + gene_len;
    print scaffold, new_start, new_end
}' > genes_2_plot.txt

cp $SCRIPT/Plot_chromosome.R ./
Rscript ./Plot_chromosome.R ${1}
rm Plot_chromosome.R
echo GOOD GENE PLOTTED

# Extract genes that blast properly to the FAD or FAD database
module load Biopython/1.83-foss-2023a
awk 'NR==FNR {ids[$1]; next} /^>/ {header=$0; id=substr($1,2); keep=(id in ids)} keep' matches.txt  unique_species_name_cds_"$FILE1" > checked_unique_species_name_cds_"$FILE1"
cat checked_unique_species_name_cds_"$FILE1"|python3 $SCRIPT/filter_too_short_long_sequence.py -o aligned_"$FILE1"_filtered -s kept_sequences -r 0.60 -m 1.2

#python3 $SCRIPT/extract_fasta_by_header.py checked_unique_species_name_cds_"$FILE1" --file kept_sequences >  tmp
mv aligned_"$FILE1"_filtered checked_unique_species_name_cds_"$FILE1"
echo LIST OF GENE BLASTING TO THE DATABASE RETRIEVED

# align sequences with muscle
module load MUSCLE/3.8.1551-GCC-9.3.0
muscle -in checked_unique_species_name_cds_"$FILE1" -out aligned_unique_species_name_cds_"$FILE1"
module load Biopython/1.83-foss-2023a
#cat aligned_unique_species_name_cds_"$FILE1"|python3 $SCRIPT/filter_too_short_long_sequence_excluding_gaps.py -o aligned_"$FILE1"_filtered -s kept_sequences -r 0.85 -m 1.15
#python3 $SCRIPT/extract_fasta_by_header.py aligned_"$FILE1"_filtered --file kept_sequences > tmp
#mv tmp aligned_"$FILE1"_filtered
echo MUSCLE ALIGNMENT DONE

# reextract the remaining unique sequences

# plot alignment
module purge
module load R/4.2.1-foss-2022a
export PATH=~/local/bin:$PATH
Rscript $SCRIPT/msa_plot.R aligned_unique_species_name_cds_"$FILE1" alignment_msa_${1}_${2}.html
echo MSA PLOT BEFORE TRIMAL GENERATED

# Trim sequence with trimal
module purge
module load trimAl/1.4.1-GCC-9.3.0

trimal -in aligned_unique_species_name_cds_"$FILE1" -gt 0.7 -cons 20 -out trimal  # -gt 0.7 more than 1-0.7=0.3 so removes sites with more than 30% of gaps. -cons 30  
#trimal -in aligned_unique_species_name_cds_"$FILE1" -gappyout -out trimal
cat  trimal | awk -v species="_${1}" '/^>/ {split($1, a, " "); print a[1] species; next} {print}' > trimal_species_name_"$FILE1"
echo TRIMAL DONE


module purge
module load R/4.2.1-foss-2022a
export PATH=~/local/bin:$PATH
Rscript $SCRIPT/msa_plot.R trimal_species_name_"$FILE1" trimal_alignment_msa_${1}_${2}.html
echo MSA PLOT AFTER TRIMAL GENERATED

# Keep unique sequences
module purge
module load Biopython/1.83-foss-2023a
python3 $SCRIPT/unique_seq.py trimal_species_name_"$FILE1" > unique_"$FILE1"
cat unique_"$FILE1" | awk -v species="_${1}" '/^>/ {split($1, a, " "); print a[1] species; next} {print}' > trimal_species_name_"$FILE1"
rm unique_"$FILE1" 

# Infer tree
#module load IQ-TREE/2.3.6-gompi-2023a
#iqtree2 -s trimal_"$FILE1" -m MFP -B 1000 --prefix $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation
	
# Plot tree with ggtree
#Rscript ./plot_tree.R --tree $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation.treefile --prefix $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation

# Clean folder
#rm matches.txt align* cds*
