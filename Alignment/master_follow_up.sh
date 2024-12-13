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
    		((counter++)) # Counting every line in the Candidate_orthogs* file
			mkdir -p $RESULTS/${PHE}/${PHE}_${counter}
			RES_PATH=$RESULTS/${PHE}/${PHE}_${counter}

			# Extract sequences for each group of paralog
    		python3 ./extract_fasta_by_header.py $RESULTS/${PHE}/All_${PHE}_combined_aa.txt --list $line | python3 ./filter_short_sequences_below_pct_mean.py -o $RES_PATH/fasta_${PHE}_"$counter"_aa -s kept_sequences -r 0.95
    		python3 ./extract_fasta_by_header.py $RESULTS/${PHE}/All_${PHE}_combined_larger_500.txt --list $line > tmp
			python3 ./extract_fasta_by_header.py tmp --file kept_sequences > $RES_PATH/fasta_${PHE}_${counter}_nt
			echo SEQUENCES FOR ${PHE} PARALOG GROUP ${counter} READY TO BE ALIGNED

			# Align, revert alignment back to nt and trim poorly aligning regions
			muscle -in $RES_PATH/fasta_${PHE}_${counter}_aa -out $RES_PATH/fasta_${PHE}_${counter}_aa.aln
    		pal2nal.pl $RES_PATH/fasta_${PHE}_${counter}_aa.aln $RES_PATH/fasta_${PHE}_${counter}_nt -output fasta -codontable 1 > $RES_PATH/fasta_${PHE}_${counter}_nt.aln
    		trimal -in $RES_PATH/fasta_${PHE}_${counter}_nt.aln -out $RES_PATH/fasta_${PHE}_${counter}_nt_trimmed.aln -automated1
			echo ALIGNMENT FOR ${PHE} ${counter} PERFORMED

			# Built ML tree and plot
			sbatch ./iqtree.sh $RES_PATH ${PHE} ${counter}
			echo SEPARATED TREE JOBS FOR ${PHE} ${counter} SUBMITTED
	done
	rm tmp
	rm kept_sequences
    echo ALL ALIGNMENTS FOR ${PHE} FINISHED
}

# Function to make a single tree for all FAD or FAR genes
Combine_all_genes_in_one_tree(){
	local PHE=$1

	# Built ML tree and plot all ${PHE} combined in a single file (so a in single tree)
	running_jobs_trees=$(squeue|grep $(whoami)| grep -P "iqtree"| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
	
	# If there are running jobs, submit with the dependency
	if [ -n "$running_jobs_trees" ]; then
    	sbatch --job-name=All${PHE} --dependency=aftercorr:$running_jobs_trees ./Combiner_all_seq_by_pheromone_family.sh $RESULTS/${PHE} ${PHE}
	else
       # No running jobs, submit without dependency
    	sbatch --job-name=All${PHE} ./Combiner_all_seq_by_pheromone_family.sh $RESULTS/${PHE} ${PHE}
	fi
}

# Main
for i in FAR; do
    echo $i
    check_pairwise_aln $i
    Ok_alignments_extract $i
    Align_candidate_orthologs $i
    Combine_all_genes_in_one_tree $i
done
