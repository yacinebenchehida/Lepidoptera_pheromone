################################
# Load libraries and set paths #
################################
module load Biopython/1.81-foss-2022b
module load BLAST+/2.14.0-gompi-2022b
module load R/4.2.1-foss-2022a
EDIRECT="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/edirect"
RESULTS="/mnt/scratch/projects/biol-specgen-2018/yacine/Pheromones/Blast/Results/"
PHEROMONE=$@

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
	mv *txt $RESULTS/$i
	rm $RESULTS/$i/db.fasta
	echo FINAL DATABASE CREATED
done

################################################################
# Create a blast data base with Heliconius melpomenes proteins #
################################################################
for i in $PHEROMONE; do
	makeblastdb -in $RESULTS/"$i"/unique_db_"$i".fasta -dbtype prot -input_type fasta -out $RESULTS/"$i"/"$i"_db -title "$i"_db
	echo BLAST DB CREATED
done

####################################################################
# Create chunks of 200 fasta sequences (to parallelize the blastx) #
####################################################################
for i in $PHEROMONE; do
	cp $RESULTS/Heliconius_melpomene_proteins.fa $RESULTS/$i
	python3 divide_fasta.py $RESULTS/$i/Heliconius_melpomene_proteins.fa
	echo PROTEINS FASTA SLICED
done

#################################################
# Run blastx on each sequence (and each chunks) #
#################################################
for i in $PHEROMONE; do
		array_number=$(ls  $RESULTS/$i/output_chunk_*|wc -l)
		sbatch --job-name=Bl"$i" --array=1-$array_number ./blastx.sh $i
	done
echo BLASTx RAN

#######################################################
# Get size of each scaffold (used later for plotting) #
#######################################################
python3 scaffold_size.py ../Inputs/Heliconius_melpomene_melpomene_Hmel2.5.scaffolds.fa |awk '$3 > 1000000' > ../Inputs/scaffold_size_information.txt

#########################################
# Combine results and extract best hits #
#########################################
running_jobs1=$(squeue|grep ybc502| grep BlFAR| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
running_jobs2=$(squeue|grep ybc502| grep BlFAD| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')

echo $(eval echo "$running_jobs1")

if [[  "$PHEROMONE" =~ "FAR" ]]; then
               sbatch --job-name=COFAR --dependency=aftercorr:$running_jobs1 ./combine_clean.sh FAR  350
                echo DATA FAR COMBINED
fi
if [[  "$PHEROMONE" =~ "FAD" ]]; then
                sbatch --job-name=COFAD --dependency=aftercorr:$running_jobs2 ./combine_clean.sh FAD 150
                echo DATA FAD COMBINED
fi

########
# Plot #
########
running_jobs3=$(squeue|grep ybc502| grep CO| awk '{print $1}'|perl -pe 's/\n/,/g'|sed 's/,$//g')
sbatch --job-name=Plot  --dependency=aftercorr:$running_jobs3 ./plot.sh 
