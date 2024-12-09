#!/bin/bash

#SBATCH --mem=5GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --time=0-00:30:00

################################
# Load libraries and set paths #
################################
module load Biopython/1.81-foss-2022b
module load BLAST+/2.14.0-gompi-2022b
module load R/4.2.1-foss-2022a
SPECIES=$@
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Alignment/Results"
INPUTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results"

# Function to create combined file for specified gene type
combine_and_align_genes() {
    local PHE=$1
    local combined_file="$RESULTS/${PHE}/All_${PHE}_combined.txt"

    # Check if the file exists and remove it if needed
    mkdir -p "$RESULTS/$PHE"
    if [ -f "$combined_file" ]; then
        rm "$combined_file"
    fi

    # Create an empty combined file
    touch "$combined_file"

    # Extract results for all species at the specified gene
    for sp in $SPECIES; do
        cat "$INPUTS/${sp}/${PHE}/${sp}_${PHE}_proteins.fasta" | awk -v species="_${sp}" '/^>/ {split($1, a, " "); print a[1] species; next} {print}' >> "$combined_file"
    done

    # Remove sequences shorter than 500 bp
    python3 ./filter_fasta_by_length.py "$combined_file" 500 > "$RESULTS/${PHE}/All_${PHE}_combined_larger_500.txt"

    # Translate nucleotide sequences to protein (for protein alignments with muscle later)
    python3 ./translate.py "$RESULTS/${PHE}/All_${PHE}_combined_larger_500.txt" > "$RESULTS/${PHE}/All_${PHE}_combined_aa.txt"

    # Alignement every pair of sequences. The script muscle.py submits a new job for each 500 pairwise comparisons.
    python3 ./muscle.py "$RESULTS/${PHE}/All_${PHE}_combined_aa.txt" $RESULTS/${PHE}/Alignments 500
}

for i in FAD; do
    combine_and_align_genes "$i"
    running_jobs_alignments=$(squeue|grep $(whoami)| grep -P "muscle_j"| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
    sbatch --job-name="$i"clean --dependency=aftercorr:$running_jobs_alignments ./run_calculate_identity.sh $RESULTS/${i}/Alignments
done
