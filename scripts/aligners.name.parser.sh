#!/bin/bash

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ usage ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# translates single letter aligner arguments into full aligner name

# ~~~~~~~~~~ synopsis ~~~~~~~~~~
# input: sequence of letters
# output: corresponding aligners

LIST=$1 # list of aligners

echo $LIST | grep -o . - | while read letter; do
  if [ "$letter" == "b" ]
  then
	echo 'bwa'
  fi
  if [ "$letter" == "s" ]
  then
	echo 'soap3-dp'
  fi
  if [ "$letter" == "n" ]
  then
	echo 'novoalign'
  fi
  if [ "$letter" == "g" ]
  then
	echo 'ngm'
  fi
  if [ "$letter" == "o" ]
  then
	echo 'bowtie2'
  fi
  if [ "$letter" == "m" ]
  then
	echo 'smalt'
  fi
done
