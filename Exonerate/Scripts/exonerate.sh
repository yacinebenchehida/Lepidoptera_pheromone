#!/bin/bash
# Author: Yacine Ben Chehida

#SBATCH --time=0-2:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=10G
#SBATCH --account=BIOL-SPECGEN-2018
#SBATCH --job-name=Exo


########################################
# 0 - Load needed modules from the hpc #
########################################
module load parallel/20230722-GCCcore-12.3.0
module load Exonerate/2.4.0-GCC-11.3.0
module load Biopython/1.83-foss-2023a
module load AUGUSTUS/3.5.0-foss-2022a

#################################################
# 1 - Set reference genomes and important paths #
#################################################
DATA="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Inputs/$1"
FASTA_REF="${1}_genome.fa"
REF_GENOME="$DATA/$FASTA_REF"
echo "Reference genome: $REF_GENOME"
export REF_GENOME  # Make visible to parallel

SPECIES=$1
export SPECIES

PHE=$2
CHUNKS_PATH="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Input/${PHE}"
RESULT="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Results/${SPECIES}/${PHE}"
RESULTS_CHUNK="${RESULT}/output_chunk_${SLURM_ARRAY_TASK_ID}"
SCRIPT="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Scripts"
MPE="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Exonerate/Scripts/extrinsic.MPE.cfg"
export SCRIPT

mkdir -p "$RESULTS_CHUNK"
cp "$CHUNKS_PATH/output_chunk_${SLURM_ARRAY_TASK_ID}.fasta" "$RESULTS_CHUNK"

#################################################
# 2 - Split each fasta entry in a separate file #
#################################################
cd "$RESULTS_CHUNK"
csplit -z -f seq_ "output_chunk_${SLURM_ARRAY_TASK_ID}.fasta" '/^>/' '{*}'
for f in seq_*; do
    name=$(grep '^>' "$f" | sed 's/^>//;s/ .*//')
    mv "$f" "${name}.fa"
done

#########################################
# 3 - Function for parallel execution   #
#########################################
run_exonerate() {
    local i="$1"
    echo "Processing FASTA: $i"
    exonerate --model protein2genome "$i" "$REF_GENOME" \
        --querytype protein --bestn 10 --showvulgar no --showtargetgff yes \
        --showalignment no --percent 20 > "${SPECIES}_${i}_exonerate_output.txt"
# extract CDSs
#    python3 "$SCRIPT/exonerate2fasta.py" \
#        -f "${SPECIES}_${i}_exonerate_output.txt" \
#        -g "$REF_GENOME" \
#        -o "${SPECIES}_${i}_tmp_cds"

# extract only genes (for augustus later)
    python3 "$SCRIPT/gene_exonerate2fasta.py" \
       -f "${SPECIES}_${i}_exonerate_output.txt" \
        -g "$REF_GENOME" \
        -o "${SPECIES}_${i}_tmp_genes"

# Get CDS directly from exonerate
    python3 "$SCRIPT/exonerate2fasta_cds.py" \
        -f "${SPECIES}_${i}_exonerate_output.txt" \
        -g "$REF_GENOME" \
        -o "${SPECIES}_${i}_cds"

# Augustus
 #   cp $MPE ./
  #  augustus --progress=true --strand=both --species=heliconius_melpomene1 "${SPECIES}_${i}_tmp_genes" --hintsfile="${SPECIES}_${i}_hints" --extrinsicCfgFile=extrinsic.MPE.cfg > ${SPECIES}_${i}_augustus_output.gff


    

}
export -f run_exonerate
# export -f run_exonerate makes the shell function available to GNU parallel,
# which runs each command in a separate subshell. Without this, parallel won't
# recognize the function name when trying to execute it.

#########################################
# 4 - Launch all jobs with parallel     #
#########################################
parallel -j 8 run_exonerate ::: *.fa

#########################################
# 5 - Merge results                     #
#########################################
#FILE1="${SPECIES}_${PHE}_combined_cds.fasta"
#FILE2="${SPECIES}_${PHE}_combined_cds.aa"

FILE1="${SPECIES}_${PHE}_combined_genes_${SLURM_ARRAY_TASK_ID}.fasta"

if [ -f "$FILE1" ]; then
    rm "$FILE1"
fi
touch "$FILE1"

#if [ -f "$FILE2" ]; then
#    rm "$FILE2"
#fi
shopt -s nullglob
files=(${SPECIES}_*_tmp_genes)

if [ ${#files[@]} -gt 0 ]; then
    cat "${files[@]}" >> "$FILE1"
    mv "$FILE1" ../
else
    echo -n "no tmp_gene files"
fi
#cat ${SPECIES}_*_tmp_cds.aa >> "$FILE2" 


# Optional: cleanup
# rm *.fa ${SPECIES}_*_exonerate_output.txt ${SPECIES}_*_tmp_cds

