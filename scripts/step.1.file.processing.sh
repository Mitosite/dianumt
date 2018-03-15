#!/bin/bash

steponestarttime=${SECONDS}

# SYNOPSIS
# input: .sam file (e.g. 1.novoalign.sam) stored in $FILENAME
# task: separates unmapped from mapped reads, separates nuclear from mitochondrial mapped reads, separates numts from non numts reads from nuclear mapped
# outputs: numts and unmapped reads as .fastq files, mitomapped and contaminating reads as sam files

FILENAME=$1
ALIGNER1=$2

echo "~~~~~~ STEP 1 ($ALIGNER1) - started processing step 1 files" >> $LOG

# extract unmapped reads into a sam file; convert to BAM then to fastq for later realignment
samtools view -hb -f4 $FILENAME | samtools bam2fq - > 1.$ALIGNER1.unmapped.fastq 2>> $LOG
echo 'step 1: extracted unmapped reads into a fastq file' >> $LOG

# extract mapped reads into a headerless sam file
samtools view -hF4 $FILENAME > 1.$ALIGNER1.mapped.sam 2>> $LOG
echo 'step 1: extracted mapped reads' >> $LOG

# extract mitochondria mapped reads into a sam file and convert it to fastq for later realignment
cat 1.$ALIGNER1.header.sam > 1.$ALIGNER1.mitomapped.sam
samtools view 1.$ALIGNER1.mapped.sam | awk '{if($3 == "MT"){print $0}}' >> 1.$ALIGNER1.mitomapped.sam 2>> $LOG
samtools view -bh 1.$ALIGNER1.mitomapped.sam | samtools bam2fq - > 1.$ALIGNER1.mitomapped.fastq  2>> $LOG
echo 'step 1: extracted mitochondrial mapped reads' >> $LOG

# extract nucleus mapped reads into a sam file
cat 1.$ALIGNER1.header.sam > 1.$ALIGNER1.nucmapped.sam
samtools view 1.$ALIGNER1.mapped.sam | awk '{if($3 != "MT"){print $0}}' >> 1.$ALIGNER1.nucmapped.sam
echo 'step 1: extracted nuclear mapped reads' >> $LOG

# sam to sorted bam to sorted bed of all nucleus mapped reads
samtools view -bh 1.$ALIGNER1.nucmapped.sam | samtools sort -o 1.$ALIGNER1.sorted.nucmapped.bam -
bedtools bamtobed -i 1.$ALIGNER1.sorted.nucmapped.bam > 1.$ALIGNER1.sorted.nucmapped.bed
rm 1.$ALIGNER1.sorted.nucmapped.bam

# extract nucleus mapped reads that intersect with numts regions (minimum of 51% of the read overlapping with a numt)
bedtools intersect -a 1.$ALIGNER1.sorted.nucmapped.bed -b $NUCLEARNUMTS -wa -u -f 0.51 > 1.$ALIGNER1.numts.bed

# extract and sort a dictionary of reads names intersecting with numt regions
cat 1.$ALIGNER1.numts.bed | cut -f4 | sort > 1.$ALIGNER1.numts.names

# remove numts from nuclear mapped bed file to obtain dictionary file of contaminating reads names
cat 1.$ALIGNER1.sorted.nucmapped.bed | cut -f4 | sort | comm -23 - 1.$ALIGNER1.numts.names > 1.$ALIGNER1.contamination.names

# re build sam file of contaminating reads using the dictionary of contaminating reads names
cat 1.$ALIGNER1.header.sam > 1.$ALIGNER1.contamination.sam
cat 1.$ALIGNER1.contamination.names | gawk -F/ '{print $1}' | fgrep -f - 1.$ALIGNER1.nucmapped.sam >> 1.$ALIGNER1.contamination.sam

# re build sam file of numts from nuclear mapped reads and convert result to fastq for later realignment
cat 1.$ALIGNER1.header.sam > 1.$ALIGNER1.numts.sam
cat 1.$ALIGNER1.numts.names | gawk -F/ '{print $1}' | fgrep -f - 1.$ALIGNER1.nucmapped.sam >> 1.$ALIGNER1.numts.sam
samtools view -bh 1.$ALIGNER1.numts.sam | samtools bam2fq - > 1.$ALIGNER1.numts.fastq 2>> $LOG
echo 'step 1: extracted numts and contaminating reads' >> $LOG

# merge relevant fastq files for steps 2 and 3 realignments
cat 1.$ALIGNER1.mitomapped.fastq 1.$ALIGNER1.numts.fastq > 1.$ALIGNER1.aligntomito.fastq 2>> $LOG
cat 1.$ALIGNER1.mitomapped.fastq 1.$ALIGNER1.numts.fastq 1.$ALIGNER1.unmapped.fastq > 1.$ALIGNER1.aligntoshiftedmito.fastq 2>> $LOG

echo 'finished processing step 1 files' >> $LOG

steponeendtimetime=${SECONDS}
steponetimedifference=`expr ${steponeendtimetime} - ${steponestarttime}`
echo "step.1.file.processing.sh; took $steponetimedifference seconds to run" >> $TIMELOG
