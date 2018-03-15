#!/bin/bash

masterstarttime=${SECONDS}

# setting up log file
export LOG=log.txt
echo 'running master shell' > $LOG

# setting up time log file
export TIMELOG=timelog.txt
echo 'elapsed time by main scripts' > $TIMELOG

# setting programs paths
export CUTADAPT=/project/soft/linux64/src/hannah/cutadapt-1.9.1-py2/bin/cutadapt # path to adapter trimming soft
export PICARD=/project/soft/linux64/src/hannah/picard/2.10.2/picard.jar # path to picard
export GATK=/project/soft/linux64/src/hannah/gatk/3.7-0/GenomeAnalysisTK.jar # path to GATK
export SAMSTAT=/project/soft/linux64/src/derekproject17/samstat/bin/samstat # path to samstat for sam files summary
export TRIMMOMATIC=/project/soft/linux64/src/derekproject17/trimmomatic/bin/trimmomatic.jar # path to trimmomatic

# setting files paths
echo 'setting paths' >> $LOG
export SHSCRIPTS=../../shscripts # shell scripts directory
export RSCRIPTS=../../Rscripts # R scripts directory
export PYSCRIPTS=../../pyscripts # python scripts directory
export DATA=../../data # path to directory with static data
export REFHUMAN=$DATA/Homo_sapiens.GRCh37.dna.primary_assembly.fa # reference human genome (hg19)
export REFMITO=$DATA/chrMT.fasta # reference mitochondrial genome (Cambridge Reference Sequence)
export REFSHIFTEDMITO=$DATA/shifted_chrMT.fasta # shifted reference mitochondrial genome (shifted by 8000 bases)
export KNOWNSNPS=$DATA/rCRS_known_sorted.vcf # reference known SNPs
export NUCLEARNUMTS=$DATA/Calabrese.et.al.2012.nuclear.bed # reference numts (nuclear positions)
export MITONUMTS=$DATA/Calabrese.et.al.2012.mitochondrial.bed # reference numts (corresponding mitochondrial position)
export INDICES=$DATA/indices # aligner indices directory

# loading required modules
source ./$SHSCRIPTS/modules.sh

# preprocessing reads: removing duplicates & trimming adapters
./$SHSCRIPTS/preprocessing.sh

######################### ALIGNING #############################

## retrieving aligners choices
export LIST1=$(cat s1_aln_choices.txt)
export LIST2=$(cat s2_aln_choices.txt)

## STEP 1: aligning input reads, GATK prep, and file processing
./$SHSCRIPTS/alignment.sh -$LIST1 -N 1 -F 0.preprocessed.reads.fastq -R human -1
wait
for ALIGNER1 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST1); do
  ./$SHSCRIPTS/quality.adjustment.sh $ALIGNER1 # indel realignment, BQSR
  ./$SHSCRIPTS/step.1.file.processing.sh 1.$ALIGNER1.preped.sam $ALIGNER1 # processing files
done
echo 'completed alignment step 1' >> $LOG

## STEP 2: remapping mitomapped and numts from step 1 to normal mitochondria
for ALIGNER1 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST1); do
  ./$SHSCRIPTS/alignment.sh -$LIST2 -N 2.$ALIGNER1 -F 1.$ALIGNER1.aligntomito.fastq -R mito -2
done
wait
echo 'completed alignment step 2' >> $LOG

## STEP 3: remapping mitomapped, numts and unmapped from step 1 to shifted mitochondria
for ALIGNER1 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST1); do
  ./$SHSCRIPTS/alignment.sh -$LIST2 -N 3.$ALIGNER1 -F 1.$ALIGNER1.aligntoshiftedmito.fastq -R shiftedmito -3
done
wait
echo 'completed alignment step 3' >> $LOG

######################### ANALYSIS #############################

## FINAL PROCESSING: variant calling, plotting & diagnostics
./$SHSCRIPTS/final.file.processing.sh
echo 'completed final file processing' >> $LOG

python $PYSCRIPTS/emailalert.py
echo 'user sent email alert' >> $LOG

masterendtimetime=${SECONDS}
mastertimedifference=`expr ${masterendtimetime} - ${masterstarttime}`
echo "master.sh: took $mastertimedifference seconds to run" >> $TIMELOG
