#!/usr/bin/env python

import pandas as pd
import sys

if len(sys.argv) < 3:
	print ("Missing files")
	exit (1)

infile1 = sys.argv[1] # 4to12
infile2 = sys.argv[2] # 12to4
outfile = sys.argv[3] # variants

# positions 4 to 12
htable = pd.read_table(infile1, low_memory= False)
wanted = htable.loc[(htable['Pos'] <= 12000)]
another_wanted = wanted.loc[(htable['Pos'] >= 4000)]

# second input file
secondtable = pd.read_table(infile2, low_memory= False)
secondtrim = secondtable.loc[(secondtable['Pos'] < 4000)] # secondtrim = 0 to 3999
thirdtrim = secondtable.loc[(secondtable['Pos'] > 12000)] # thirdtrim = 12000-16669
stopnow = thirdtrim.loc[(thirdtrim['Pos'] <= 16569)] # stopnow = 12000-16569 - trimming non existent bases

# concatenate into 1 to 16569
concat = pd.concat([secondtrim, another_wanted])
secondconcat = pd.concat([concat, stopnow])

# arrange index to start from 1 to make equal to the positions
secondconcat.index = pd.RangeIndex(start=1, stop=16570, step=1)
secondconcat.index.name = "Position" # naming index column

# retrieving meaningful columns from concatenated mitocaller output
secondconcat.drop(secondconcat.columns[0:2],axis=1, inplace=True)
secondconcat.drop(secondconcat.columns[1:9],axis=1, inplace=True)
secondconcat.drop(secondconcat.columns[2:13],axis=1, inplace=True)
secondconcat.drop(secondconcat.columns[2:9],axis=1, inplace=True)
secondconcat.drop(secondconcat.columns[5:],axis=1, inplace=True)

# removing column headers present at each line
secondconcat['FilteredDepth'] = secondconcat['FilteredDepth'].map(lambda x: x.lstrip('DepthFilter:'))
secondconcat['Ref'] = secondconcat['Ref'].map(lambda x: x.lstrip('RefBase:'))
secondconcat['Genotype'] = secondconcat['Genotype'].map(lambda x: x.lstrip('Genotype'))
secondconcat['Genotype'] = secondconcat['Genotype'].map(lambda x: x.lstrip(':'))
secondconcat['AlleleFractions'] = secondconcat['AlleleFractions'].map(lambda x: x.lstrip('Frequency:'))

# extract non matches between reference and called base at the position
secondconcat.loc[(secondconcat['Ref'] != secondconcat['Genotype'])]

# sorting by reference allele less than 95%
newtable = secondconcat.loc[(secondconcat['MajorAlleleFraction'] < 0.9500)]

# export variants to csv
newtable.to_csv(outfile)
