#!/usr/bin/env python

import pandas as pd
import sys

if len(sys.argv) < 1:
  print("Missing mitocaller output file.")
  exit(1)

# input file
infile = sys.argv[1]
table = pd.read_table(infile, low_memory= False)

# retrieving meaningful columns from concatenated mitocaller output
table.drop(table.columns[0:2],axis=1, inplace=True)
table.drop(table.columns[1:9],axis=1, inplace=True)
table.drop(table.columns[2:13],axis=1, inplace=True)
table.drop(table.columns[2:9],axis=1, inplace=True)
table.drop(table.columns[5:],axis=1, inplace=True)
table.index = pd.RangeIndex(start=1, stop=16570, step=1)
table.index.name = "Position"

# removing column headers present at each line
table['FilteredDepth'] = table['FilteredDepth'].map(lambda x: x.lstrip('DepthFilter:'))
table['Ref'] = table['Ref'].map(lambda x: x.lstrip('RefBase:'))
table['Genotype'] = table['Genotype'].map(lambda x: x.lstrip('Genotype'))
table['Genotype'] = table['Genotype'].map(lambda x: x.lstrip(':'))
table['AlleleFractions'] = table['AlleleFractions'].map(lambda x: x.lstrip('Frequency:'))

# extract non matches between reference and called base at the position
table.loc[(table['Ref'] != table['Genotype'])]

# sorting by reference allele less than 95%
newtable = table.loc[(table['MajorAlleleFraction'] < 0.9500)]

# export variants to csv
name = str(infile) + ".variants.summary.csv"
newtable.to_csv(name)
