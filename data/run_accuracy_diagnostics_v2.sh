##this a shell script to run the accuracy and diagnostic .sh for many different files
# to be run in the same directory as your files and the accuracy and diagnostics scripts

##load modules if needed
#module load bedtools
#module load derek
#module load samtools
#module load qualimap


BOWTIE= $1
BWA=$2
NGM=$3
NOVOALIGN=$4
SMALT=$5

for VAR in 0.1 1.0 2.0; do
	./accuracy.sh $VAR$BOWTIE $VAR_bowtie_accuracy.txt 10
	./diagnostics.sh $VAR$BOWTIE $VAR_bowtie_diagnostics.txt 
done

for VAR in 0.1 1.0 2.0; do
	./accuracy.sh $VAR$BWA $VAR_bwa_accuracy.txt 10
	./diagnostics.sh $VAR$BWA $VAR_bwa_diagnostics.txt 
done

for VAR in 0.1 1.0 2.0; do
	./accuracy.sh $VAR$NGM $VAR_ngm_accuracy.txt 10
	./diagnostics.sh $VAR$NGM $VAR_ngm_diagnostics.txt 
done

for VAR in 0.1 1.0 2.0; do
	./accuracy.sh $VAR$NOVOALIGN $VAR_novoalign_accuracy.txt 10
	./diagnostics.sh $VAR$NOVOALIGN $VAR_novoalign_diagnostics.txt 
done

for VAR in 0.1 1.0 2.0; do
	./accuracy.sh $VAR$SMALT $VAR_smalt_accuracy.txt 10
	./diagnostics.sh $VAR$SMALT $VAR_smalt_diagnostics.txt
done

