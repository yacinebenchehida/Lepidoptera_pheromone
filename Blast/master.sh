################################
# Load libraries and set paths #
################################
module load Biopython/1.81-foss-2022b
module load BLAST+/2.14.0-gompi-2022b
EDIRECT="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/edirect"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results"

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
touch $RESULTS/db_FAR.fasta

for i in uniprotkb uniparc uniref
do
	curl -L "https://rest.uniprot.org/${i}/stream?download=true&format=fasta&includeIsoform=true&query=%28%22fatty+acyl-coa+reductase%22+AND+%28taxonomy_id%3A7088%29%29" >> $RESULTS/db_FAR.fasta
done
echo UNIPROT DATA DOWNLOADED

################################
# download data from gene bank #
################################
$EDIRECT/esearch -db Protein -query "fatty acyl-coa reductase[All Fields] AND Lepidoptera[Organism]"| $EDIRECT/efetch -format fasta >> $RESULTS/db_FAR.fasta
echo GENBANK DATA DOWNLOADED

#####################################################################
# Keep unique entries in the constitute database (removed doublons) #
#####################################################################
python3 unique_fasta.py $RESULTS/db_FAR.fasta FAR > $RESULTS/unique_db_FAR.fasta
mv *txt $RESULTS
rm $RESULTS/db_FAR.fasta

echo FINAL DATABASE CREATED

################################################################
# Create a blast data base with Heliconius melpomenes proteins #
################################################################
makeblastdb -in $RESULTS/unique_db_FAR.fasta -dbtype prot -input_type fasta -out $RESULTS/FAR_db -title FAR_db
echo BLAST DB CREATED

####################################################################
# Create chunks of 200 fasta sequences (to parallelize the blastp) #
####################################################################
python3 divide_fasta.py ../Results/Heliconius_melpomene_proteins.fa
echo HELICONIUS PROTEINS FASTA SLICED

#################################################
# Run blastp on each sequence (and each chunks) #
#################################################
array_number=$(ls  $RESULTS/output_chunk_*|wc -l)

for CHUNK in $RESULTS/chunk*
do
	sbatch --array=1-$array_number ./blastp.sh
done

echo BLASTP RAN
