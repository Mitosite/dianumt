#!/bin/bash

pairedpreprostarttime=${SECONDS}

echo 'started input fastq file preprocessing' >> $LOG

# retrieving adapter sequence
ADAPTER=adapters.fa

# trimming adapters
echo 'running trimmomatic' >> $LOG
java -jar $TRIMMOMATIC PE -phred33 input1.fastq input2.fastq paired_trimmed_1.fastq.gz unpaired_trimmed_1.fastq.gz paired_trimmed_2.fastq.gz unpaired_trimmed_2.fastq.gz ILLUMINACLIP:$ADAPTER:2:2:2 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:15 MINLEN:20 2>> $LOG
echo 'trimmed with trimmomatic' >> $LOG

# removing duplicates
echo 'running superdeduper' >> $LOG
super_deduper -1 paired_trimmed_1.fastq.gz -2 paired_trimmed_2.fastq.gz -p paired_deduped_trimmed 2>> $LOG
mv paired_deduped_trimmed_nodup_PE1.fastq 0.preprocessed.readsfile1.fastq
mv paired_deduped_trimmed_nodup_PE2.fastq 0.preprocessed.readsfile2.fastq
echo 'removed duplicates with superdeduper' >> $LOG

echo 'completed preprocessing' >> $LOG

pairedpreproendtimetime=${SECONDS}
pairedpreprotimedifference=`expr ${pairedpreproendtimetime} - ${pairedpreprostarttime}`
echo "preprocessing.pairedend.sh; took $pairedpreprotimedifference seconds to run" >> $TIMELOG
