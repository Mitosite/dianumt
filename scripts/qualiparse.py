# parses Qualimap txt output for easy access to key information

import xlwt
import getopt, sys, os


# INITIALISATION

def usage():
	print ('''\n  Qualimap txt file parsing ~~\n
  Required options:
    -i		input file
    -o		output file name with .xls extension\n
    ''')

try:
  opts, args = getopt.getopt(sys.argv[1:], 'i:o:')
except getopt.GetoptError as err:
  print(str('  ' + str(err)))
  sys.exit()

if len(opts) != 2:
  usage()
  sys.exit()

for o,a in opts:
  if o == '-i': input_file = str(a)
  if o == '-o': output_file = str(a)
  else: assert 'unhandled option'


# MAIN INSTRUCTIONS

# open file and read as list
with open(input_file, 'r') as in_file:
    data = in_file.readlines()

# set up the strings to search
searches = ['number of reads',
           'number of mapped reads',
           'number of duplicated reads',
           'duplication rate',
           'general error rate',
           'number of mismatches',
           'number of insertions',
           'number of bases',
           'number of deletions']

# will store the content of the line for each successful search
lines = ['']*len(searches)

# extract relevant line for each of the searches string
for i in range (0, len(searches)):
    lines[i] = [line.rstrip() for line in data if searches[i] in line]

# avoid error if could not find search string
for i in range (0, len(searches)):
    if lines[i] == []:
        lines[i] = ['no string matching this search; N/A']

# fill in relevant variables
reads = lines[0][0].split(' ')[len(lines[0][0].split(' '))-1]
reads = int(reads.replace(',',''))
mapped = lines[1][0].split(' ')[len(lines[1][0].split(' '))-2]
mapped = int(mapped.replace(',',''))
unmapped = reads - mapped
estimated_duplicate = lines[2][0].split(' ')[len(lines[2][0].split(' '))-1]
estimated_duplicate = int(estimated_duplicate.replace(',',''))
duplication_rate = lines[3][0].split(' ')[len(lines[3][0].split(' '))-1]
error_rate = lines[4][0].split(' ')[len(lines[4][0].split(' '))-1]
mismatches = lines[5][0].split(' ')[len(lines[5][0].split(' '))-1]
mismatches = mismatches.replace(',','')
insertions = lines[6][0].split(' ')[len(lines[6][0].split(' '))-1]
insertions = insertions.replace(',','')
bases = lines[7][0].split(' ')[len(lines[7][0].split(' '))-2]
bases = int(bases.replace(',',''))
deletions = lines[8][0].split(' ')[len(lines[8][0].split(' '))-1]
deletions = deletions.replace(',','')

# extract mapped chromosomes data
for j in range(0, len(data)):
    if data[j].rstrip() == '>>>>>>> Coverage per contig':
        chrs = data[j+2:j+27]

# set up sorted chromosome order
chromosomes = []
chromosomes.extend(range(1,23))
chromosomes.extend(['X','Y','MT'] )

length = [] # total number of bases for each chromosome
mapped_bases = [] # number of bases mapped to each chromosome
for k in range(0, len(chromosomes)):
    for chromosome in chrs:
        if str(chromosome.split('\t')[1]) == str(chromosomes[k]):
            length.extend([chromosome.split('\t')[2]])
            mapped_bases.extend([chromosome.split('\t')[3]])

# extract a list of chromosomes to which reads have been mapped
hit_chromosomes = []
for l in range(0, len(chromosomes)):
    if int(mapped_bases[l]) != 0:
        hit_chromosomes.extend([chromosomes[l]])

# exporting in Excel spreadsheet
wb = xlwt.Workbook()
ws = wb.add_sheet('main')

# ws.write(line, column, item), with line and column starting from 0

ws.write(1, 1, output_file)

ws.write(3, 1, 'total bases (bp)'); ws.write(3, 2, str(bases))
ws.write(4, 1, 'total reads n°'); ws.write(4, 2, str(reads))

ws.write(6, 1, 'mapped reads n°'); ws.write(6, 2, str(mapped)); ws.write(6, 3, str(str(round((mapped/reads)*100,2)) + '%'))
ws.write(7, 1, 'unmapped reads n°'); ws.write(7, 2, str(unmapped)); ws.write(7, 3, str(str(round((unmapped/reads)*100,2)) + '%'))
ws.write(8, 1, 'hit chromosomes'); ws.write(8, 2, str(hit_chromosomes))

ws.write(10, 1,'error rate'); ws.write(10, 2, str(error_rate))
ws.write(11, 1,'estimated duplicated reads'); ws.write(11, 2, str(estimated_duplicate))
ws.write(12, 1,'duplication rate'); ws.write(12, 2, str(duplication_rate))
ws.write(13, 1,'mismatches'); ws.write(13, 2, str(mismatches))
ws.write(14, 1,'insertions'); ws.write(14, 2, str(insertions))
ws.write(15, 1,'deletions'); ws.write(15, 2, str(deletions))

ws.write(17, 1,'detailed mapping per chromosome')
ws.write(18, 1,'name'); ws.write(18, 2, 'length'); ws.write(18, 3, 'number of mapped bases');

count = 0
for m in range(0, len(chromosomes)):
    if int(mapped_bases[m]) != 0:
        ws.write(int(19+count),1,str(chromosomes[m]))
        ws.write(int(19+count),2,str(length[m]))
        ws.write(int(19+count),3,str(mapped_bases[m]))
        count+=1

wb.save(output_file)

# end
