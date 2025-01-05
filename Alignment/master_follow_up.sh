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
	python3 ./filter_unique_species.py $RESULTS/${PHE}/Candidate_orthogs_${PHE}.txt $RESULTS/${PHE}/Candidate_orthogs_${PHE}_filtered_uniq_sp.txt
    echo CANDIDATE ORTHOLOGS FILE GENERATED
}

# Function to make a fasta file for all the genes aligning, revert alignment back to nucleotids, trim excessive gaps, and realign
Align_candidate_orthologs(){
	module purge
	module load PAL2NAL/14-GCCcore-10.3.0
	module load MUSCLE/3.8.1551-GCC-9.3.0
	module load trimAl/1.4.1-GCC-9.3.0
	module load R/4.2.1-foss-2022a
	module load Biopython/1.81-foss-2022b
	module load seqtk/1.3-GCC-11.3.0

	local PHE=$1
	export PATH=~/local/bin:$PATH # To use pandoc (plotting msa)

	cat $RESULTS/${PHE}/Candidate_orthogs_${PHE}_filtered_uniq_sp.txt |while read line; do \
    		((counter++)) # Counting every line in the Candidate_orthogs* file
			mkdir -p $RESULTS/${PHE}/${PHE}_${counter}
			RES_PATH=$RESULTS/${PHE}/${PHE}_${counter}

			# Extract sequences for each group of paralog
    		python3 ./extract_fasta_by_header.py $RESULTS/${PHE}/All_${PHE}_combined_aa.txt --list $line | python3 ./filter_too_short_long_sequence.py -o $RES_PATH/fasta_${PHE}_"$counter"_aa -s kept_sequences -r 0.95 -m 1.2
    		python3 ./extract_fasta_by_header.py $RESULTS/${PHE}/All_${PHE}_combined_larger_500.txt --list $line > tmp
			python3 ./extract_fasta_by_header.py tmp --file kept_sequences > $RES_PATH/fasta_${PHE}_${counter}_nt
			echo SEQUENCES FOR ${PHE} PARALOG GROUP ${counter} READY TO BE ALIGNED

			# Align, revert alignment back to nt and trim poorly aligning regions
			muscle -in $RES_PATH/fasta_${PHE}_${counter}_aa -out $RES_PATH/fasta_${PHE}_${counter}_aa.aln
			perl -pe 's/-/./g' $RES_PATH/fasta_${PHE}_${counter}_aa.aln > cleaned_alignment_${counter}.fasta
			Rscript ./msa_plot.R cleaned_alignment_${counter}.fasta $RES_PATH/${PHE}_${counter}_aa_alignment.html
			rm cleaned_alignment_${counter}.fasta
    		pal2nal.pl $RES_PATH/fasta_${PHE}_${counter}_aa.aln $RES_PATH/fasta_${PHE}_${counter}_nt -output fasta -codontable 1 > $RES_PATH/fasta_${PHE}_${counter}_nt.aln
    		trimal -in $RES_PATH/fasta_${PHE}_${counter}_nt.aln -out $RES_PATH/fasta_${PHE}_${counter}_nt_trimmed.aln -automated1
			python3 ./Pick_unique_seq_per_sp.py $RES_PATH/fasta_${PHE}_${counter}_nt_trimmed.aln $RES_PATH/fasta_${PHE}_${counter}_nt_trimmed_unique_seq_sp.aln
			Rscript ./msa_plot.R $RES_PATH/fasta_${PHE}_${counter}_nt_trimmed_unique_seq_sp.aln $RES_PATH/${PHE}_${counter}_nt_trimmed_alignment.html
			echo ALIGNMENT FOR ${PHE} ${counter} PERFORMED

			# Built ML tree and plot
			sbatch ./iqtree.sh $RES_PATH ${PHE} ${counter}
			echo SEPARATED TREE JOBS FOR ${PHE} ${counter} SUBMITTED
	done
	rm tmp kept_sequences
    echo ALL ALIGNMENTS FOR ${PHE} FINISHED
}

# Function to make a single tree for all FAD or FAR genes
Combine_all_genes_in_one_tree_by_adding_gaps(){
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

Combine_all_genes_in_one_tree_for_manual_curation(){
	local PHE=$1

	module purge
	module load PAL2NAL/14-GCCcore-10.3.0
	module load MUSCLE/3.8.1551-GCC-9.3.0
	module load trimAl/1.4.1-GCC-9.3.0
	module load R/4.2.1-foss-2022a
	module load Biopython/1.81-foss-2022b

	# Extract from the FAD FAR gene only the one with at least 3 taxa
	python3 ./Prepare_fasta_4_manual_curation.py $RESULTS/${PHE} ${PHE}
	# Get a combined AA sequence for all fAD/FAR genes with at least 3 taxa
	python3 ./extract_fasta_by_header.py $RESULTS/${PHE}/All_${PHE}_combined_aa.txt --file $RESULTS/${PHE}/${PHE}_manual_curation/Merged_manual_curation_${PHE}.txt > $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_aa_4_manual_curation.fasta
	# Get a combined NT sequence for all fAD/FAR genes with at least 3 taxa
	python3 ./extract_fasta_by_header.py $RESULTS/${PHE}/All_${PHE}_combined_larger_500.txt --file $RESULTS/${PHE}/${PHE}_manual_curation/Merged_manual_curation_${PHE}.txt > $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_nt_4_manual_curation.fasta

	# AA alignment with muscle
	muscle -in $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_aa_4_manual_curation.fasta -out $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_aa_4_manual_curation_aln.fasta
	# Plot AA alignment
	perl -pe 's/-/./g' $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_aa_4_manual_curation_aln.fasta > $RESULTS/${PHE}/${PHE}_manual_curation/tmp_fasta_4_plotting_msa
	export PATH=~/local/bin:$PATH # To use pandoc (plotting msa)
	Rscript ./msa_plot.R $RESULTS/${PHE}/${PHE}_manual_curation/tmp_fasta_4_plotting_msa $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_aa_4_manual_curation_aln.html
	rm $RESULTS/${PHE}/${PHE}_manual_curation/tmp_fasta_4_plotting_msa
	
	# Revert alignment back to nucleotide
	pal2nal.pl $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_aa_4_manual_curation_aln.fasta $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_nt_4_manual_curation.fasta -output fasta -codontable 1 > $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_NT_4_manual_curation_aln.fasta
	# Make sure all sequences are unique (if not pick randomly one sequence)
	python3 ./unique_sequence.py $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_NT_4_manual_curation_aln.fasta $RESULTS/${PHE}/${PHE}_manual_curation/unique_${PHE}_manual_curation_NT_4_manual_curation_aln.fasta
	# Trim alignment with trimal 
	trimal -in  $RESULTS/${PHE}/${PHE}_manual_curation/unique_${PHE}_manual_curation_NT_4_manual_curation_aln.fasta -out  $RESULTS/${PHE}/${PHE}_manual_curation/trimmed_unique_${PHE}_manual_curation_NT_4_manual_curation_aln.fasta  -gt 0.5 -cons 50
	# Plot final alignment before tree
	Rscript ./msa_plot.R $RESULTS/${PHE}/${PHE}_manual_curation/trimmed_unique_${PHE}_manual_curation_NT_4_manual_curation_aln.fasta $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation_NT_4_manual_curation_aln.html

	# ML tree with IQtree
	module load IQ-TREE/2.3.6-gompi-2023a
	iqtree2 -s $RESULTS/${PHE}/${PHE}_manual_curation/trimmed_unique_${PHE}_manual_curation_NT_4_manual_curation_aln.fasta -m MFP -B 1000 --prefix $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation
	# Plot tree with ggtree
	Rscript ./plot_tree.R --tree $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation.treefile --prefix $RESULTS/${PHE}/${PHE}_manual_curation/${PHE}_manual_curation

	# Clean folder with all the intermediate files
	rm $RESULTS/${PHE}/${PHE}_manual_curation/*fasta 

}

# Main
for i in FAD FAR; do
	echo $i
    check_pairwise_aln $i
    Ok_alignments_extract $i
    Align_candidate_orthologs $i
	Combine_all_genes_in_one_tree_by_adding_gaps $i
	Combine_all_genes_in_one_tree_for_manual_curation $i
done
