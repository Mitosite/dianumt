#!/bin/bash

qualstarttime=${SECONDS}

ALIGNER1=$1

samtools view -bh 1.$ALIGNER1.sam > 1.$ALIGNER1.bam
samtools sort -o sorted.$ALIGNER1.bam 1.$ALIGNER1.bam
samtools index sorted.$ALIGNER1.bam

# Add Read Groups 
java -jar $PICARD AddOrReplaceReadGroups I=sorted.$ALIGNER1.bam O=sorted.$ALIGNER1.read.bam RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20 2>> $LOG
echo 'added read groups' >> $LOG

#### Indel realignments
# First step - RealignerTargetCreator
samtools index sorted.$ALIGNER1.read.bam
java -jar $GATK -T RealignerTargetCreator -R $REFHUMAN -I sorted.$ALIGNER1.read.bam -known $KNOWNSNPS -U ALLOW_SEQ_DICT_INCOMPATIBILITY -o realigner.$ALIGNER1.intervals 2>> $LOG
#### Second step - IndelRealigner
java -jar $GATK -T IndelRealigner -R $REFHUMAN -I sorted.$ALIGNER1.read.bam -U ALLOW_SEQ_DICT_INCOMPATIBILITY -known $KNOWNSNPS -targetIntervals realigner.$ALIGNER1.intervals -o realigned.$ALIGNER1.bam 2>> $LOG
echo 'completed indel realignment' >> $LOG

#### Base Recalibration
# first step - BaseRecalibrator
samtools index realigned.$ALIGNER1.bam
java -jar $GATK -T BaseRecalibrator -R $REFHUMAN -I realigned.$ALIGNER1.bam -knownSites $KNOWNSNPS -o recal.$ALIGNER1.table -U ALLOW_SEQ_DICT_INCOMPATIBILITY 2>> $LOG
# second step - PrintReads
java -jar $GATK -T PrintReads -R $REFHUMAN -I realigned.$ALIGNER1.bam -BQSR recal.$ALIGNER1.table -o recal.$ALIGNER1.bam -U ALLOW_SEQ_DICT_INCOMPATIBILITY 2>> $LOG
echo 'completed base recalibration' >> $LOG

# final file
samtools view -h recal.$ALIGNER1.bam > 1.$ALIGNER1.preped.sam

qualendtimetime=${SECONDS}
qualtimedifference=`expr ${qualendtimetime} - ${qualstarttime}`
echo "quality.adjustment.sh; took $qualtimedifference seconds to run" >> $TIMELOG
