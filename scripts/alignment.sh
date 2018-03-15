#!/bin/bash

alignstarttime=${SECONDS}

# set default behaviour (run no aligner)
b=FALSE
s=FALSE
n=FALSE
g=FALSE
o=FALSE
m=FALSE

# parse arguments
while getopts ":N:F:R:bsngom123" opt; do
  case $opt in
    N) # desired output file name, to which will be added the name of the aligner and .sam
      NAME=$OPTARG
      ;;
    F) # .fastq file location / path
      FILEPATH=$OPTARG
      ;;
    R) # reference genome path
      REFERENCE=$OPTARG
      ;;
    b) # bwa
      b=TRUE
      ;;
    s) # soap3-dp
      s=TRUE
      ;;
    n) # novoalign
      n=TRUE
      ;;
    g) # ngm
      g=TRUE
      ;;
    o) # bowtie2
      o=TRUE
      ;;
    m) # smalt
      m=TRUE
      ;;
    1) # step 1
      STEP=1
      ;;
    2) # step 2
      STEP=2
      ;;
    3) # step 3
      STEP=3
      ;;
    \?) # error
      echo "invalid option: -$OPTARG" >> $LOG
      ;;
  esac
done

echo "aligning reads... (alignment step $STEP of 3)" >> $LOG

# BWA
if [ "$b" == "TRUE" ]
then
	echo "bwa is true" >> $LOG
	bwa mem $INDICES/bwa/$REFERENCE.bwa.fa $FILEPATH > $NAME.bwa.sam 2>> $LOG
	samtools view -H $NAME.bwa.sam > $STEP.bwa.header.sam
else
	echo "bwa is false" >> $LOG
fi

# SOAP3-DP
if [ "$s" == "TRUE" ]
then
	echo "soap3-dp is true" >> $LOG
	cp $INDICES/soap3-dp/*.ini .
	soap3-dp single $INDICES/soap3-dp/$REFERENCE.index $FILEPATH -L 100 -b 2 2>> $LOG
	rm *.ini
	samtools merge $NAME.soap3-dp.sam $FILEPATH.*
	rm $FILEPATH.*
	samtools view -H $NAME.soap3-dp.sam > $STEP.soap3-dp.header.sam
else
	echo "soap3-dp is false" >> $LOG
fi

# NOVOALIGN
if [ "$n" == "TRUE" ]
then
	echo "novoalign is true" >> $LOG
	novoalign -d $INDICES/novoalign/$REFERENCE.novoalign.index -o SAM -f $FILEPATH > $NAME.novoalign.sam 2>> $LOG
	samtools view -H $NAME.novoalign.sam > $STEP.novoalign.header.sam
else
	echo "novoalign is false" >> $LOG
fi

# NGM
if [ "$g" == "TRUE" ]
then
	echo "ngm is true" >> $LOG
	ngm -q $FILEPATH -r $INDICES/ngm/$REFERENCE.fasta -o $NAME.ngm.sam -t 4 2>> $LOG
	samtools view -H $NAME.ngm.sam > $STEP.ngm.header.sam
else
	echo "ngm is false" >> $LOG
fi

# BOWTIE2
if [ "$o" == "TRUE" ]
then
	echo "bowtie2 is true" >> $LOG
	bowtie2 -x $INDICES/bowtie2/$REFERENCE.bowtie2.index -U $FILEPATH -S $NAME.bowtie2.sam 2>> $LOG
	samtools view -H $NAME.bowtie2.sam > $STEP.bowtie2.header.sam
else
	echo "bowtie2 is false" >> $LOG
fi

# SMALT
if [ "$m" == "TRUE" ]
then
	echo "smalt is true" >> $LOG
	smalt map -f samsoft -o $NAME.smalt.sam $INDICES/smalt/$REFERENCE.smalt $FILEPATH 2>> $LOG
	echo "completed alignment with smalt" >> $LOG
	samtools view -H $NAME.smalt.sam > $STEP.smalt.header.sam
else
	echo "smalt is false" >> $LOG
fi

alignendtimetime=${SECONDS}
aligntimedifference=`expr ${alignendtimetime} - ${alignstarttime}`
echo "alignment.sh; took $aligntimedifference seconds to run" >> $TIMELOG
