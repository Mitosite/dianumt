#!/bin/bash

finalstarttime=${SECONDS}

# SYNOPSIS
# input: sam files for step 1 and steps 2+3
# tasks: process sam files for variant calling, diagnostics and plotting
# outputs: mitocaller outputs for before and after pipeline



# ~~~~~~~~~~~~~~~~~~~~~~~~~~ VARIANT CALLING ~~~~~~~~~~~~~~~~~~~~~~~~~~
# uses mitocaller to call variants, extracts call different from reference

echo 'calling variants on step 1 files' >> $LOG
for ALIGNER1 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST1); do
  samtools view -bh 1.$ALIGNER1.sam > 1.$ALIGNER1.bam
  samtools sort -o 1.$ALIGNER1.bam 1.$ALIGNER1.bam
  samtools index 1.$ALIGNER1.bam
  samtools view -bh 1.$ALIGNER1.bam MT > 1.$ALIGNER1.mt.bam # extract only mitochondrial reads
  echo 'base recalibration' >> $LOG
  samtools calmd -Abr 1.$ALIGNER1.mt.bam $REFMITO > 1.$ALIGNER1.baqed.bam
  echo 'running mitocaller' >> $LOG
  mitoCaller -m -b 1.$ALIGNER1.baqed.bam -r $REFMITO > 1.$ALIGNER1.mitocaller
  rm 1.$ALIGNER1.bam* 1.$ALIGNER1.mt.bam 1.$ALIGNER1.baqed.bam # remove temporary files
done

echo 'calling variants on steps 2 and 3 files' >> $LOG
for ALIGNER1 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST1); do
  for ALIGNER2 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST2); do

    # call variants on normal mitochondria genome at positions 4000 to 12000 (step 2)
    samtools view -bh 2.$ALIGNER1.$ALIGNER2.sam > 2.$ALIGNER1.$ALIGNER2.bam
    samtools sort -o 2.$ALIGNER1.$ALIGNER2.bam 2.$ALIGNER1.$ALIGNER2.bam
    samtools index 2.$ALIGNER1.$ALIGNER2.bam
    samtools view -bh 2.$ALIGNER1.$ALIGNER2.bam MT:4000-12000 > 4to12.$ALIGNER1.$ALIGNER2.bam
    echo 'base recalibration' >> $LOG
    samtools calmd -Abr 4to12.$ALIGNER1.$ALIGNER2.bam $REFMITO > baq.4to12.$ALIGNER1.$ALIGNER2.bam
    echo 'running mitocaller' >> $LOG
    mitoCaller -m -b baq.4to12.$ALIGNER1.$ALIGNER2.bam -r $REFMITO > 4to12.$ALIGNER1.$ALIGNER2.mitocaller
    rm 2.$ALIGNER1.$ALIGNER2.bam 2.$ALIGNER1.$ALIGNER2.bam.bai 4to12.$ALIGNER1.$ALIGNER2.bam baq.4to12.$ALIGNER1.$ALIGNER2.bam

    # call variants on shifted mitochondria genome at positions 1 to 3999 and 12001 to 16569 (step 3)
    samtools view -H 2.$ALIGNER1.$ALIGNER2.sam > offset.$ALIGNER1.$ALIGNER2.sam
    samtools view 3.$ALIGNER1.$ALIGNER2.sam | awk '{$4 = ((($4 + 8000)-1) % 16569)+1; print}' - | sed -e 's/ /\t/g' >> offset.$ALIGNER1.$ALIGNER2.sam
    samtools view -bh offset.$ALIGNER1.$ALIGNER2.sam > offset.$ALIGNER1.$ALIGNER2.bam
    samtools sort -o offset.$ALIGNER1.$ALIGNER2.bam offset.$ALIGNER1.$ALIGNER2.bam
    samtools index offset.$ALIGNER1.$ALIGNER2.bam
    samtools view -b offset.$ALIGNER1.$ALIGNER2.bam MT:12001-16569 > 12to0.$ALIGNER1.$ALIGNER2.bam
    samtools view -b offset.$ALIGNER1.$ALIGNER2.bam MT:0-3999 > 0to4.$ALIGNER1.$ALIGNER2.bam
    samtools merge 12to4.$ALIGNER1.$ALIGNER2.bam 12to0.$ALIGNER1.$ALIGNER2.bam 0to4.$ALIGNER1.$ALIGNER2.bam
    echo 'base recalibration' >> $LOG
    samtools calmd -Abr 12to4.$ALIGNER1.$ALIGNER2.bam $REFMITO > baq.12to4.$ALIGNER1.$ALIGNER2.bam 2> /dev/null
    echo 'running mitocaller' >> $LOG
    mitoCaller -m -b baq.12to4.$ALIGNER1.$ALIGNER2.bam -r $REFMITO > 12to4.$ALIGNER1.$ALIGNER2.mitocaller
    rm offset.$ALIGNER1.$ALIGNER2.sam offset.$ALIGNER1.$ALIGNER2.bam offset.$ALIGNER1.$ALIGNER2.bam.bai
    rm 12to0.$ALIGNER1.$ALIGNER2.bam 0to4.$ALIGNER1.$ALIGNER2.bam 12to4.$ALIGNER1.$ALIGNER2.bam baq.12to4.$ALIGNER1.$ALIGNER2.bam

    # merge mitocaller outputs: combine 4to12 and 12to4, output coverage table and variants summary
    echo 'merging steps 2 and 3 mitocaller outputs, outputting coverage and variants summaries' >> $LOG
    ./$PYSCRIPTS/merging.and.interpreting.mitocaller.py 4to12.$ALIGNER1.$ALIGNER2.mitocaller 12to4.$ALIGNER1.$ALIGNER2.mitocaller 2+3.$ALIGNER1.$ALIGNER2.mitocaller.variants.summary.csv
    ./$PYSCRIPTS/merging.and.interpreting.mitocaller.py 4to12.$ALIGNER1.$ALIGNER2.mitocaller 12to4.$ALIGNER1.$ALIGNER2.mitocaller 2+3.$ALIGNER1.$ALIGNER2.mitocaller.variants.summary.csv

  done
done
echo 'completed variant calling' >> $LOG



# ~~~~~~~~~~~~~~~~~~~~~~~~~ IDENTIFY BEST ALIGNER COMBINATION ~~~~~~~~~~~~~~~~~~~~~~~~~~
# identify best combination based on criterion of maximum mitochondria coverage

export COVSUMMARY=coverage.summary.txt
echo -e 'total depth\taligner1\taligner2' > $COVSUMMARY

echo 'starting steps 2 and 3 files processing for identifying best aligner combination' >> $LOG
for ALIGNER1 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST1); do
  for ALIGNER2 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST2); do

    # extract total mitochondria coverage atfer pipeline and export to file
    ./$SHSCRIPTS/merge.coverage.sh $ALIGNER1 $ALIGNER2
    TOTALCOV=$(cat 2+3.$ALIGNER1.$ALIGNER2.coverage | awk 'BEGIN{}{sum+=$3} END{print sum}')
    echo -e "$TOTALCOV\t$ALIGNER1\t$ALIGNER2" >> $COVSUMMARY

  done
done

echo 'made mitochondria coverage file from mitocaller output of steps 2+3' >> $LOG
export BESTALIGNER1=$(sed '1d' coverage.summary.txt | sort -k 1 | head -n1 | awk '{{print $2}}')
export BESTALIGNER2=$(sed '1d' coverage.summary.txt | sort -k 1 | head -n1 | awk '{{print $3}}')
echo "the best results for this job were found by running $BESTALIGNER1 at step 1 and $BESTALIGNER2 at step 2" >> $LOG



# ~~~~~~~~~~~~~~~~~~~~~~~~~ FILE PROCESSING ~~~~~~~~~~~~~~~~~~~~~~~~~~
# makes files required for diagnostics and plotting numts and contamination data

echo 'starting step 1 files processing for analysis and plotting' >> $LOG

# running samstat
$SAMSTAT 1.$BESTALIGNER1.sam
echo 'ran samstat on step 1 file' >> $LOG

# extract per-base filtered depth for step 1 file
cat 1.$BESTALIGNER1.mitocaller | sed '1d' | awk {'print $1"\t"$2"\t"$4'} | sed 's/Depth://' > 1.$BESTALIGNER1.coverage
echo 'made mitochondria coverage file from mitocaller output of step 1' >> $LOG

# make bed file of contamination (sam to sorted bam to bed)
samtools view -bh 1.$BESTALIGNER1.contamination.sam | samtools sort -o 1.$BESTALIGNER1.sorted.contamination.bam -
bedtools bamtobed -i 1.$BESTALIGNER1.sorted.contamination.bam > 1.$BESTALIGNER1.contamination.bed
rm 1.$BESTALIGNER1.sorted.contamination.bam

# extract unique names of numts reads were mapped to
bedtools intersect -a 1.$BESTALIGNER1.sorted.nucmapped.bed -b $NUCLEARNUMTS -wb -f 0.51 | cut -f10 | sort | uniq > 1.$BESTALIGNER1.foundnumts.names

# make bed file of numts as found in step 1 (nuclear coordinates for plotting)
fgrep -f 1.$BESTALIGNER1.foundnumts.names $NUCLEARNUMTS > 1.$BESTALIGNER1.nuclearnumts.bed

# make bed file of numts as found in step 1 (mitochondrial coordinates for depth analysis)
fgrep -f 1.$BESTALIGNER1.foundnumts.names $MITONUMTS > 1.$BESTALIGNER1.mitonumts.bed

# ~~~~ prepare relevant files for plotting ~~~~

# prepare mitochondrial positions numts bed file for plotting (add * as directionality)
mv 1.$BESTALIGNER1.mitonumts.bed temporarymitonumts.bed
awk '$5=="" {$5="*"}1' OFS="\t" temporarymitonumts.bed > 1.$BESTALIGNER1.mitonumts.bed
rm temporarymitonumts.bed

# prepare numts bed file for plotting
mv 1.$BESTALIGNER1.numts.bed temporarynumts.bed
cut -d$'\t' -f1-4,6 temporarynumts.bed > 1.$BESTALIGNER1.numts.bed # delete 5th column containing score
rm temporarynumts.bed

# prepare contamination bed file for plotting
mv 1.$BESTALIGNER1.contamination.bed temporarycontamination.bed
cut -d$'\t' -f1-4,6 temporarycontamination.bed > 1.$BESTALIGNER1.contamination.bed # delete 5th column containing score
rm temporarycontamination.bed

# compile bed file of all reads to plot with relevant keyword in the name column, and chr for start of chromosome
cat 1.$BESTALIGNER1.numts.bed | awk '{print "chr"$1"\t"$2"\t"$3"\tnumtsreads\t*"}' > $BESTALIGNER1.readsforplotting.bed
cat 1.$BESTALIGNER1.mitonumts.bed | awk '{print "chrM\t"$2"\t"$3"\tmitonumts\t"$5}' >> $BESTALIGNER1.readsforplotting.bed
cat 1.$BESTALIGNER1.nuclearnumts.bed | awk '{print "chr"$1"\t"$2"\t"$3"\tnuclearnumts\t"$5}' >> $BESTALIGNER1.readsforplotting.bed
cat 1.$BESTALIGNER1.contamination.bed | awk '{print "chr"$1"\t"$2"\t"$3"\tcontamination\t"$5}' >> $BESTALIGNER1.readsforplotting.bed

# delete GL sequences from file
mv $BESTALIGNER1.readsforplotting.bed temporaryreads.bed
cat temporaryreads.bed | sed '/^chrGL/d' > $BESTALIGNER1.readsforplotting.bed
rm temporaryreads.bed

echo 'completed file processing for best aligner combination' >> $LOG



# ~~~~~~~~~~~~~~~~~~~~~~~~~ DATA ANALYSIS ~~~~~~~~~~~~~~~~~~~~~~~~~~
# collects data about the best aligner outputs and exports it to summary file

# create summary file for the current combination of step 1 and step 2 aligners
export SUMMARY=summary.txt

# writing job details to summary output
echo -e "~~ Job details ~~" > $SUMMARY
INPUTREADS=$(bc <<< "scale=0; $(cat input1.fastq | wc -l)/4")
OUTPUTREADS=$(bc <<< "scale=0; $(cat 0.preprocessed.reads.fastq | wc -l)/4")
echo -e "Read count in input file: $INPUTREADS"  >> $SUMMARY
if [ -s adapter.txt ]
then
    ADAPTER=$(cat adapter.txt)
    echo -e "Single end reads pipeline - adapter used: $ADAPTER"  >> $SUMMARY
else
    ADAPTERS=$(cat adapters.fa)
    echo -e "Paired end reads pipeline at step 1 ~~\nAdapters: $ADAPTERS"  >> $SUMMARY
fi
echo -e "Read count in preprocessed file: $OUTPUTREADS" >> $SUMMARY
STEP1ALIGNERS=$(echo $(./$SHSCRIPTS/aligners.name.parser.sh $LIST1) | sed  's/\n/ /' | sed 's/ /, /')
STEP2ALIGNERS=$(echo $(./$SHSCRIPTS/aligners.name.parser.sh $LIST2) | sed  's/\n/ /' | sed 's/ /, /')
echo -e "Aligner(s) used in step 1: $STEP1ALIGNERS" >> $SUMMARY
echo -e "Aligner(s) used in step 2: $STEP2ALIGNERS" >> $SUMMARY
echo -e "The best results for this job were found by running $BESTALIGNER1 at step 1 and $BESTALIGNER2 at step 2" >> $SUMMARY
wait

# STEP 1
echo -e "\n~~ Step 1 (align to whole human genome) ~~" >> $SUMMARY
CONTAMINATION=$(samtools view -c 1.$BESTALIGNER1.contamination.sam)
echo -e "Estimated number of contaminating reads: $CONTAMINATION\n" >> $SUMMARY
TOTAL=$(samtools view -c 1.$BESTALIGNER1.sam)
MAPPED=$(samtools view -c 1.$BESTALIGNER1.mapped.sam)
MITOMAPPED=$(samtools view -c 1.$BESTALIGNER1.mitomapped.sam)
NUMTSMAPPED=$(samtools view -c 1.$BESTALIGNER1.numts.sam)
UNMAPPED=$(bc <<< "scale=0; $(cat 1.$BESTALIGNER1.unmapped.fastq | wc -l)/4")
echo -e "Total reads: $TOTAL\nTotal mapped reads: $MAPPED\nMapped to mitochondrion: $MITOMAPPED\nMapped to numts: $NUMTSMAPPED\nUnmapped: $UNMAPPED" >> $SUMMARY
UNIQUENUMTS=$(cat 1.$BESTALIGNER1.foundnumts.names | wc -l)
echo -e "Number of unique numts reads were mapped to: $UNIQUENUMTS" >> $SUMMARY
SUM1=$(cat 1.$BESTALIGNER1.coverage | awk '{sum+=$3;} END {print sum}')
AVG1=$(cat 1.$BESTALIGNER1.coverage | awk '{sum+=$3;} END {print sum/NR}')
STDEV1=$(cat 1.$BESTALIGNER1.coverage | awk '{sum+=$3; sumsq+=$3*$3} END {print sqrt(sumsq/NR - (sum/NR)**2)}')
echo -e "\nDepth of coverage in the entire mitochondrion:\nSum: $SUM1\nAverage: $AVG1\nStandard deviation: $STDEV1" >> $SUMMARY
# isolate total, mean, and standard deviation of depth at locations of numts found in step 1, in step 1 file
MITOSUM1=$(./$PYSCRIPTS/depth.at.numts.py 1.$BESTALIGNER1.mitonumts.bed 1.$BESTALIGNER1.coverage | awk '{sum+=$2;} END {print sum}')
MITOAVG1=$(./$PYSCRIPTS/depth.at.numts.py 1.$BESTALIGNER1.mitonumts.bed 1.$BESTALIGNER1.coverage | awk '{range+=$1; sum+=$2;} END {print sum/range}')
MITOSTDEV1=$(./$PYSCRIPTS/depth.at.numts.py 1.$BESTALIGNER1.mitonumts.bed 1.$BESTALIGNER1.coverage | awk '{range+=$1; sum+=$2; sumsq+=$2*$2} END {print sqrt(sumsq/range - (sum/range)**2)}')
echo -e "\nDepth of coverage in the mitochondrion at numts:\nSum: $MITOSUM1\nAverage: $MITOAVG1\nStandard deviation: $MITOSTDEV1" >> $SUMMARY
wait

# STEP 2
echo -e "\n~~ Step 2 (align to mitochondria only) ~~" >> $SUMMARY
SUM2=$(cat 2+3.$BESTALIGNER1.$BESTALIGNER2.coverage | awk '{sum+=$3;} END {print sum}')
AVG2=$(cat 2+3.$BESTALIGNER1.$BESTALIGNER2.coverage | awk '{sum+=$3;} END {print sum/NR}')
STDEV2=$(cat 2+3.$BESTALIGNER1.$BESTALIGNER2.coverage | awk '{sum+=$3; sumsq+=$3*$3} END {print sqrt(sumsq/NR - (sum/NR)**2)}')
wait
echo -e "\nDepth of coverage in the entire mitochondrion:\nSum: $SUM2\nAverage: $AVG2\nStandard deviation: $STDEV2" >> $SUMMARY
# isolate total, mean, and standard deviation of depth at locations of numts found in step 1, in steps 2+3 file
MITOSUM2=$(./$PYSCRIPTS/depth.at.numts.py 1.$BESTALIGNER1.mitonumts.bed 2+3.$BESTALIGNER1.$BESTALIGNER2.coverage | awk '{sum+=$2;} END {print sum}')
MITOAVG2=$(./$PYSCRIPTS/depth.at.numts.py 1.$BESTALIGNER1.mitonumts.bed 2+3.$BESTALIGNER1.$BESTALIGNER2.coverage | awk '{range+=$1; sum+=$2;} END {print sum/range}')
MITOSTDEV2=$(./$PYSCRIPTS/depth.at.numts.py 1.$BESTALIGNER1.mitonumts.bed 2+3.$BESTALIGNER1.$BESTALIGNER2.coverage | awk '{range+=$1; sum+=$2; sumsq+=$2*$2} END {print sqrt(sumsq/range - (sum/range)**2)}')
echo -e "\nDepth of coverage in the mitochondrion at numts:\nSum: $MITOSUM2\nAverage: $MITOAVG2\nStandard deviation: $MITOSTDEV2" >> $SUMMARY
wait

# extract variants from the mitocaller output
./$PYSCRIPTS/interpreting.mitocaller.py 1.$BESTALIGNER1.mitocaller
./$PYSCRIPTS/interpreting.mitocaller.py 1.$BESTALIGNER1.mitocaller
echo 'exported summary tables of variants for steps 1 and 2' >> $LOG

# integrate R scripts for plotting data
Rscript $RSCRIPTS/mitograph.R $BESTALIGNER1.readsforplotting.bed 1.$BESTALIGNER1.coverage 2+3.$BESTALIGNER1.$BESTALIGNER2.coverage
echo 'plotted karyotypes' >> $LOG

echo 'completed data analysis and plotting' >> $LOG



# ~~~~~~~~~~~~~~~~~~~~~~~~~ DOWNLOADS LIST ~~~~~~~~~~~~~~~~~~~~~~~~~~
# create a zipped file to be downloaded from the HAMSTR results webpage for the job

python -c "import sys; sys.path.insert(0, '../../pyscripts'); from mitositefuncs import zipper; zipper('output.zip', 'summary.txt', 'log.txt', 'coverage.summary.txt', 's1_aln_choices.txt', 's2_aln_choices.txt', '1.$BESTALIGNER1.sam.samstat.html', '1.$BESTALIGNER1.sam', '2.$BESTALIGNER1.$BESTALIGNER2.sam', '3.$BESTALIGNER1.$BESTALIGNER2.sam', '1.$BESTALIGNER1.aligntomito.fastq', '1.$BESTALIGNER1.aligntoshiftedmito.fastq', '1.$BESTALIGNER1.mitocaller', '4to12.$BESTALIGNER1.$BESTALIGNER2.mitocaller', '12to4.$BESTALIGNER1.$BESTALIGNER2.mitocaller', '1.$BESTALIGNER1.numts.sam', '0.preprocessed.reads.fastq', '1.$BESTALIGNER1.mitocaller.variants.summary.csv', '2+3.$BESTALIGNER1.$BESTALIGNER2.mitocaller.variants.summary.csv')"
echo 'created zipped downloads file' >> $LOG

echo 'completed final file processing' >> $LOG

finalendtimetime=${SECONDS}
finaltimedifference=`expr ${finalendtimetime} - ${finalstarttime}`
echo "final.file.processing.sh; took $finaltimedifference seconds to run" >> $TIMELOG
