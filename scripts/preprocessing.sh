#!/bin/bash

preprostarttime=${SECONDS}

echo 'started input fastq file preprocessing' >> $LOG

# retrieving adapter sequence
ADAPTER=$(cat adapter.txt)

# trimming adapters
echo 'running curadapt' >> $LOG
$CUTADAPT -a $ADAPTER -q 25 input1.fastq > trimmed.input1.fastq 2>> $LOG
echo 'trimmed adapters with cutadapt' >> $LOG

# removing duplicates
echo 'running superdeduper' >> $LOG
super_deduper -s 1 -l 50 -U trimmed.input1.fastq -p preprocessed_reads 2>> /dev/null
mv preprocessed_reads_nodup_PE1.fastq 0.preprocessed.reads.fastq
echo 'removed duplicates with superdeduper' >> $LOG

# removing temporary file
rm trimmed.input1.fastq

echo 'completed preprocessing' >> $LOG

preproendtimetime=${SECONDS}
preprotimedifference=`expr ${preproendtimetime} - ${preprostarttime}`
echo "preprocessing.sh; took $preprotimedifference seconds to run" >> $TIMELOG
