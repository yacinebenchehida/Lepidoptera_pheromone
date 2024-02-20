################################
# Load libraries and set paths #
################################
module load Biopython/1.81-foss-2022b
EDIRECT="/mnt/scratch/projects/biol-specgen-2018/yacine/Tools/edirect"

###############################
# Download CDS from melpomene #
###############################
wget http://download.lepbase.org/v4/sequence/Heliconius_melpomene_melpomene_Hmel2.5.cds.fa.gz

##############################
# download data from uniprot #
##############################
touch db_FAR.fasta

for i in uniprotkb uniparc uniref
do
	curl -L "https://rest.uniprot.org/${i}/stream?download=true&format=fasta&includeIsoform=true&query=%28%22fatty+acyl-coa+reductase%22+AND+%28taxonomy_id%3A7088%29%29" >> db_FAR.fasta
done


################################
# download data from gene bank #
################################
$EDIRECT/esearch -db Protein -query "fatty acyl-coa reductase[All Fields] AND Lepidoptera[Organism]"| $EDIRECT/efetch -format fasta >> db_FAR.fasta

#####################################################################
# Keep unique entries in the constitute database (removed doublons) #
#####################################################################
python3 unique_fasta.py db_FAR.fasta FAR > unique_db_FAR.fasta
