#!/bin/bash

mergestarttime=${SECONDS}

ALIGNER1=$1
ALIGNER2=$2

# extract coverage data
cat 4to12.$ALIGNER1.$ALIGNER2.mitocaller | awk '{print $1"\t"$2"\t"$4}' | sed 's/Depth://' > 4to12.$ALIGNER1.$ALIGNER2.coverage
cat 12to4.$ALIGNER1.$ALIGNER2.mitocaller | awk '{print $1"\t"$2"\t"$4}' | sed 's/Depth://' > 12to4.$ALIGNER1.$ALIGNER2.coverage

# extract positions 1 to 3999
for POS in {1..3999}
do
   awk -v x=$POS '{if ($2 == x){print $0}}' 12to4.$ALIGNER1.$ALIGNER2.coverage >> 2+3.$ALIGNER1.$ALIGNER2.coverage
done

# extract positions 4000 to 12000
for POS in {4000..12000}
do
   awk -v x=$POS '{if ($2 == x){print $0}}' 4to12.$ALIGNER1.$ALIGNER2.coverage >> 2+3.$ALIGNER1.$ALIGNER2.coverage
done

# extract positions 12001 to 16569
for POS in {12001..16569}
do
   awk -v x=$POS '{if ($2 == x){print $0}}' 12to4.$ALIGNER1.$ALIGNER2.coverage >> 2+3.$ALIGNER1.$ALIGNER2.coverage
done

# remove files
rm 4to12.$ALIGNER1.$ALIGNER2.coverage 12to4.$ALIGNER1.$ALIGNER2.coverage

mergeendtimetime=${SECONDS}
mergetimedifference=`expr ${mergeendtimetime} - ${mergestarttime}`
echo "merge.coverage.sh; took $mergetimedifference seconds to run" >> $TIMELOG
