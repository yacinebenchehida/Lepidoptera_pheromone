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
EDIRECT="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/edirect"
PHEROMONE=$(echo $@ | awk '{print $1"\t"$2}')
SPECIES=$(echo $@ | awk '{print $3}')
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results/$SPECIES"
INPUTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Inputs/$SPECIES"

################################
# Download CDS in right format #
################################
mkdir -p $RESULTS
cat $INPUTS/*cds.fa| perl -ne 's/[^\x00-\x7F]+/ /g; print;'|perl -pe 's/\.1//g' > $RESULTS/"$SPECIES"_proteins.fa


##############################
# download data from uniprot #
##############################
if [[ "$PHEROMONE" =~ "FAR" ]]; then
    mkdir -p  $RESULTS/FAR
    touch $RESULTS/FAR/db.fasta	
    for i in uniprotkb uniparc uniref
	do
  		curl -L "https://rest.uniprot.org/${i}/stream?download=true&format=fasta&includeIsoform=true&query=%28%22fatty+acyl-coa+reductase%22+AND+%28taxonomy_id%3A7088%29%29" >> $RESULTS/FAR/db.fasta
		echo "UNIPROT DATA FROM $i DOWNLOADED"
	done
fi

if [[ "$PHEROMONE" =~ "FAD" ]]; then
    mkdir -p $RESULTS/FAD
    touch $RESULTS/FAD/db.fasta
    for i in uniprotkb uniparc uniref
	do
		curl -L "https://rest.uniprot.org/${i}/stream?download=true&format=fasta&includeIsoform=true&query=%28%22fatty+acyl-CoA+desaturase%22+AND+%28taxonomy_id%3A7088%29%29" >> $RESULTS/FAD/db.fasta
		echo "UNIPROT DATA FROM $i DOWNLOADED"
	done
fi

################################
# download data from gene bank #
################################
if [[ "$PHEROMONE" =~ "FAR" ]]; then
		$EDIRECT/esearch -db Protein -query "fatty acyl-coa reductase[All Fields] AND Lepidoptera[Organism]"| $EDIRECT/efetch -format fasta >> $RESULTS/FAR/db.fasta
		echo GENBANK DATA DOWNLOADED
fi
if [[ "$PHEROMONE" =~ "FAD" ]]; then
		$EDIRECT/esearch -db Protein -query "fatty acyl-coa desaturase[All Fields] AND Lepidoptera[Organism]"| $EDIRECT/efetch -format fasta >> $RESULTS/FAD/db.fasta
		echo GENBANK DATA DOWNLOADED
fi

#####################################################################
# Keep unique entries in the constitute database (removed doublons) #
#####################################################################
for i in $PHEROMONE; do
	python3 unique_fasta.py $RESULTS/$i/db.fasta "$i" > $RESULTS/$i/unique_db_"$i".fasta
	rm $RESULTS/$i/db.fasta
	echo FINAL DATABASE CREATED
done

##########################################
# Create a blast data base with proteins #
##########################################
for i in $PHEROMONE; do
	makeblastdb -in $RESULTS/"$i"/unique_db_"$i".fasta -dbtype prot -input_type fasta -out $RESULTS/"$i"/"$i"_db -title "$i"_db
	echo BLAST DB CREATED
done

####################################################################
# Create chunks of 200 fasta sequences (to parallelize the blastx) #
####################################################################
for i in $PHEROMONE; do
	cp $RESULTS/"$SPECIES"_proteins.fa $RESULTS/$i
	python3 divide_fasta.py $RESULTS/$i/"$SPECIES"_proteins.fa
	echo PROTEINS FASTA SLICED
done

#################################################
# Run blastx on each sequence (and each chunks) #
#################################################
for i in $PHEROMONE; do
		array_number=$(ls  $RESULTS/$i/output_chunk_*|wc -l)
		SP=$(echo $SPECIES|cut -c 1-3)
		NAME=$(echo "$SP""$i")
		sbatch --job-name="$NAME" --array=1-$array_number ./blastx.sh $SPECIES $i
	done
echo BLASTx RAN

#######################################################
# Get size of each scaffold (used later for plotting) #
#######################################################
python3 scaffold_size.py ../Inputs/$SPECIES/*genome.fa |awk '$3 > 1000000' > ../Inputs/$SPECIES/scaffold_size_information.txt

#########################################
# Combine results and extract best hits #
#########################################
SP=$(echo $SPECIES|cut -c 1-3)
echo -e "${SP}FAR" 
echo -e "${SP}FAD"

running_jobs1=$(squeue|grep ybc502| grep -P "${SP}FAR"| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
running_jobs2=$(squeue|grep ybc502| grep -P "${SP}FAD"| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')

echo $(eval echo "$running_jobs1")

if [[  "$PHEROMONE" =~ "FAR" ]]; then
               sbatch --job-name="$SP"_R --dependency=aftercorr:$running_jobs1 ./combine_clean.sh $SPECIES FAR  500
                echo DATA FAR COMBINED
fi
if [[  "$PHEROMONE" =~ "FAD" ]]; then
                sbatch --job-name="$SP"_D --dependency=aftercorr:$running_jobs2 ./combine_clean.sh $SPECIES FAD 250
                echo DATA FAD COMBINED
fi

########
# Plot #
########
running_jobs3=$(squeue|grep ybc502| grep -P "${SP}_" | awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
echo $running_jobs3
sbatch --job-name=Plot  --dependency=aftercorr:$running_jobs3 ./plot.sh $SPECIES
