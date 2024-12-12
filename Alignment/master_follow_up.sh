#!/bin/bash

#SBATCH --mem=5GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-04:00:00

RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Alignment/Results"
INPUTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results"

# Function to check pairwise alignment quality
check_pairwise_aln() {
	module load Biopython/1.81-foss-2022b
	local PHE=$1
	python3 ./calculate_identity.py $RESULTS/${PHE}/Alignments
	echo IDENTITY CALCULATED
}

# Function to extract set of genes that aligned correctly
Ok_alignments_extract() {
	module purge
	module load Cython/3.0.8-GCCcore-12.2.0
	local PHE=$1
	python -c "import filter_orthologs; filter_orthologs.main()" $RESULTS/${PHE}/Alignments 300 0.70 3 > $RESULTS/${PHE}/Candidate_orthogs_${PHE}.txt
    echo CANDIDATE ORTHOLOGS FILE GENERATED
}

# Function to make a fasta file for all the genes aligning, revert alignment back to nucleotids, trim excessive gaps, and realign
Align_candidate_orthologs(){
	module purge
	module load PAL2NAL/14-GCCcore-10.3.0
	module load MUSCLE/3.8.1551-GCC-9.3.0
	module load trimAl/1.4.1-GCC-9.3.0
	
	local PHE=$1
	cat $RESULTS/${PHE}/Candidate_orthogs_${PHE}.txt |while read line; do \
    		((counter++))
    		python3 ./extract_fasta_by_header.py $RESULTS/${PHE}/All_${PHE}_combined_aa.txt --list $line | python3 ./filter_short_sequences_below_pct_mean.py -o $RESULTS/${PHE}/fasta_${PHE}_"$counter"_aa -s kept_sequences -r 0.95
    		python3 ./extract_fasta_by_header.py $RESULTS/${PHE}/All_${PHE}_combined_larger_500.txt --list $line > tmp
			python3 ./extract_fasta_by_header.py tmp --file kept_sequences > $RESULTS/${PHE}/fasta_${PHE}_${counter}_nt
			muscle -in $RESULTS/${PHE}/fasta_${PHE}_${counter}_aa -out $RESULTS/${PHE}/fasta_${PHE}_${counter}_aa.aln
    		pal2nal.pl $RESULTS/${PHE}/fasta_${PHE}_${counter}_aa.aln $RESULTS/${PHE}/fasta_${PHE}_${counter}_nt -output fasta -codontable 1 > $RESULTS/${PHE}/fasta_${PHE}_${counter}_nt.aln
    		trimal -in $RESULTS/${PHE}/fasta_${PHE}_${counter}_nt.aln -out $RESULTS/${PHE}/fasta_${PHE}_${counter}_nt_trimmed.aln -automated1
	done
    echo ALIGNMENT FOR ${PHE} FINISHED
}

# Main
for i in FAD FAD; do
	echo $i
    check_pairwise_aln $i
    Ok_alignments_extract $i
    Align_candidate_orthologs $i
done
