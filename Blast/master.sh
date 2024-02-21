################################
# Load libraries and set paths #
################################
module load Biopython/1.81-foss-2022b
module load BLAST+/2.14.0-gompi-2022b
EDIRECT="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/edirect"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results/$1"
PHEROMONE=$1

mkdir -p $RESULTS

###############################
# Download CDS from melpomene #
###############################
wget http://download.lepbase.org/v4/sequence/Heliconius_melpomene_melpomene_Hmel2.5.cds.fa.gz
mv Heliconius_melpomene_melpomene_Hmel2.5.cds.fa.gz $RESULTS
echo MELPOMENE PROTEINS DOWNLOADED
zcat $RESULTS/Heliconius_melpomene_melpomene_Hmel2.5.cds.fa.gz| perl -ne 's/[^\x00-\x7F]+/ /g; print;' > $RESULTS/Heliconius_melpomene_proteins.fa
rm $RESULTS/Heliconius_melpomene_melpomene_Hmel2.5.cds.fa.gz

##############################
# download data from uniprot #
##############################
touch $RESULTS/db_"$PHEROMONE".fasta

if [ "$PHEROMONE" = "FAR" ]; then
    for i in uniprotkb uniparc uniref
	do
  		curl -L "https://rest.uniprot.org/${i}/stream?download=true&format=fasta&includeIsoform=true&query=%28%22fatty+acyl-coa+reductase%22+AND+%28taxonomy_id%3A7088%29%29" >> $RESULTS/db_"$PHEROMONE".fasta
		echo "UNIPROT DATA FROM $i DOWNLOADED"
	done
elif [ "$PHEROMONE" = "FAD" ]; then
    for i in uniprotkb uniparc uniref
	do
		curl -L "https://rest.uniprot.org/${i}/stream?download=true&format=fasta&includeIsoform=true&query=%28%22fatty+acyl-CoA+desaturase%22+AND+%28taxonomy_id%3A7088%29%29" >> $RESULTS/db_"$PHEROMONE".fasta
		echo "UNIPROT DATA FROM $i DOWNLOADED"
	done
else
    echo "unknown pheromone"
fi

################################
# download data from gene bank #
################################
if [ "$PHEROMONE" = "FAR" ]; then
		$EDIRECT/esearch -db Protein -query "fatty acyl-coa reductase[All Fields] AND Lepidoptera[Organism]"| $EDIRECT/efetch -format fasta >> $RESULTS/db_"$PHEROMONE".fasta
		echo GENBANK DATA DOWNLOADED
elif [ "$PHEROMONE" = "FAD" ]; then
		$EDIRECT/esearch -db Protein -query "fatty acyl-coa desaturase[All Fields] AND Lepidoptera[Organism]"| $EDIRECT/efetch -format fasta >> $RESULTS/db_"$PHEROMONE".fasta
		echo GENBANK DATA DOWNLOADED
else
    echo "unknown pheromone"
fi

#####################################################################
# Keep unique entries in the constitute database (removed doublons) #
#####################################################################
python3 unique_fasta.py $RESULTS/db_"$PHEROMONE".fasta "$PHEROMONE" > $RESULTS/unique_db_"$PHEROMONE".fasta
mv *txt $RESULTS
rm $RESULTS/db_"$PHEROMONE".fasta

echo FINAL DATABASE CREATED

################################################################
# Create a blast data base with Heliconius melpomenes proteins #
################################################################
makeblastdb -in $RESULTS/unique_db_"$PHEROMONE".fasta -dbtype prot -input_type fasta -out $RESULTS/"$PHEROMONE"_db -title "$PHEROMONE"_db
echo BLAST DB CREATED

####################################################################
# Create chunks of 200 fasta sequences (to parallelize the blastx) #
####################################################################
python3 divide_fasta.py $RESULTS/Heliconius_melpomene_proteins.fa
echo HELICONIUS PROTEINS FASTA SLICED

#################################################
# Run blastx on each sequence (and each chunks) #
#################################################
array_number=$(ls  $RESULTS/output_chunk_*|wc -l)

for CHUNK in $RESULTS/chunk*
do
	sbatch --array=1-$array_number ./blastx.sh $PHEROMONE
done

echo BLASTx RAN

#########################################
# Combine results and extract best hits #
#########################################
running_jobs1=$(squeue|grep ybc502| grep blastx| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
echo $(eval echo "$running_jobs1")
sbatch --job-name=COMB --dependency=aftercorr:$running_jobs1 ./combine_clean.sh $PHEROMONE
