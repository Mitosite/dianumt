#!/bin/bash

pairedalignstarttime=${SECONDS}

# set default behaviour (run no aligner)
b=FALSE
s=FALSE
n=FALSE
g=FALSE
o=FALSE
m=FALSE

# parse arguments
while getopts ":N:F:H:R:bsngom123" opt; do
  case $opt in

    N) # input file name
      export NAME=$OPTARG
      ;;
    F) # .fastq file location / path
      export FILENAME1=$OPTARG
      ;;
    H) # need both files for paired end reads
      export FILENAME2=$OPTARG
      ;;
    R) # reference genome path
      export REFERENCE=$OPTARG
      ;;
    b) # bwa
      b=TRUE
      ;;
    s) # soap3-dp
      s=FALSE
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
      echo "step 1" >&2
      export STEP=1
      ;;
    2) # step 2
      echo "step 2" >&2
      export STEP=2
      ;;
    3) # step 3
      echo "step 3" >&2
      export STEP=3
      ;;
    \?) # error
      echo "invalid option: -$OPTARG" >> $LOG
      ;;
  esac
done

echo "aligning paired end reads... (alignment step $STEP of 3)" >> $LOG

# BWA
if [ "$b" == "TRUE" ]
then
	echo "bwa is true" >> $LOG
	bwa mem $INDICES/bwa/$REFERENCE.bwa.fa $FILENAME1 $FILENAME2 > $NAME.bwa.sam 2>> $LOG
	samtools view -H $NAME.bwa.sam > $STEP.bwa.header.sam 
else
	echo "bwa is false" >> $LOG
fi

# SOAP3-DP
if [ "$s" == "TRUE" ]
then
	echo "soap3-dp is true" >> $LOG
	cp $INDICES/soap3-dp/*.ini .
	soap3-dp pair $INDICES/soap3-dp/$REFERENCE.index $FILENAME1 $FILENAME2 -L 100 -b 2 2>> $LOG
	rm *.ini
	samtools merge $NAME.soap3-dp.sam $FILENAME1.*
	rm $FILENAME1.*
	samtools view -H $NAME.soap3-dp.sam > $STEP.soap3-dp.header.sam
else
	echo "soap3-dp is false" >> $LOG
fi

# BOWTIE2
if [ "$o" == "TRUE" ]
then
	echo "bowtie2 is true" >> $LOG
	bowtie2 -x $INDICES/bowtie2/$REFERENCE.bowtie2.index -1 $FILENAME1 -2 $FILENAME2 -S $NAME.bowtie2.sam 2>> $LOG
	samtools view -H $NAME.bowtie2.sam > $STEP.bowtie2.header.sam
else
	echo "bowtie2 is false" >> $LOG
fi

# NOVOALIGN
if [ "$n" == "TRUE" ]
then
	echo "novoalign is true" >> $LOG
	novoalign -d $INDICES/novoalign/$REFERENCE.fasta -o SAM -f $FILENAME1 $FILENAME2 > $NAME.novoalign.sam 2>> $LOG
	samtools view -H $NAME.novoalign.sam > $STEP.novoalign.header.sam
else
	echo "novoalign is false" >> $LOG
fi

# NGM
if [ "$g" == "TRUE" ]
then
	echo "ngm is true" >> $LOG
	ls $INDICES/ngm
	ngm -r $INDICES/ngm/$REFERENCE.ngm.* -1 $FILENAME1 -2 $FILENAME2 -o $NAME.ngm.sam -t 4 2>> $LOG
	samtools view -H $NAME.ngm.sam > $STEP.ngm.header.sam
else
	echo "ngm is false" >> $LOG
fi

# SMALT
if [ "$m" == "TRUE" ]
then
	smalt map -f samsoft -o $NAME.smalt.sam $INDICES/smalt/$REFERENCE.smalt $FILENAME1 $FILENAME2 2>> $LOG
	samtools view -H $NAME.smalt.sam > $STEP.smalt.header.sam
else
	echo "smalt is false" >> $LOG
fi

pairedalignendtimetime=${SECONDS}
pairedaligntimedifference=`expr ${pairedalignendtimetime} - ${pairedalignstarttime}`
echo "alignment.pairedend.sh; took $pairedaligntimedifference seconds to run" >> $TIMELOG
