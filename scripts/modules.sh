#!/bin/sh -i

module load derek/1.0
module load superdeduper
module load bwa
module load samtools
module load bedtools
module load soap3-dp
module load bowtie
module load mitocaller
module load R
module load python

echo 'modules successfully loaded' >> $LOG
