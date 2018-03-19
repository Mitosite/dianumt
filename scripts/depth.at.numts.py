#!/usr/bin/env python

import sys

if len(sys.argv) < 2:
  print("Usage : coverage  bed_file  coverage_file")
  print("E.g.  : coverage 1.smalt.mitonumts.bed 1.smalt.coverage")
  exit(1)

bed_file = sys.argv[1]
cov_file = sys.argv[2]

fb = open(bed_file, 'r')
bed_data = fb.readlines()
fb.close()

fc = open(cov_file, 'r')
cov_data = fc.readlines()
fc.close()

cov_key1 = []
cov_key2 = []
cov_key3 = []
cov_count = 0
for line in cov_data:
   line.strip()
   columns = line.split()
   cov_count = cov_count + 1 # number of lines in coverage file
   cov_key1.append(columns[0])
   cov_key2.append(columns[1])
   cov_key3.append(columns[2])

for line in bed_data:
   line.strip()
   columns    = line.split()
   line_mt    = columns[0]
   line_start = int(columns[1])
   line_end   = int(columns[2])
   line_key   = columns[3]

   sum = 0.0
   for i in range(line_start,line_end+1):
      sum = sum + int(cov_key3[i-1])
   print(str(line_end-line_start) + "\t" + str(sum)) # prints number of bases and total coverage of the region
