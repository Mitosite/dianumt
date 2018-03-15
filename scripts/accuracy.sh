#### implements aligners accuracy measures as described in: Ruffalo, Matthew, Thomas LaFramboise, and Mehmet Koyutürk. "Comparative analysis of algorithms for next-generation sequencing read alignment." Bioinformatics 27.20 (2011): 2790-2796.

# this shell must be used on sam files with namings consistent with code; i.e. reads mapped to the mitochondria must show "MT" in column 5; reads from the mitochondria must have an underscore-delimited name starting with "MT"; and reads from the nucleus must have an underscore-delimited name starting with anything but MT (e.g. origin chromosome number).

# setting paths and variables
export DATA=../../data # path to directory with static data
export NUCLEARNUMTS=$DATA/Calabrese.et.al.2012.nuclear.bed # reference numts (nuclear positions)
export SHSCRIPTS=../../shscripts # shell scripts directory
export LOG=log.txt
THRESHOLD=$1

# checking input mapq
if [ -z "$1" ]
  then
    echo -e "no mapq threshold argument provided\n"
    exit 1
fi

# loading modules
source ./$SHSCRIPTS/modules.sh

# retrieving aligners choices at first step
export LIST1=$(cat s1_aln_choices.txt)

for ALIGNER1 in $(./$SHSCRIPTS/aligners.name.parser.sh $LIST1); do

# setting input and output file names for the current iteration
FILENAME=1.$ALIGNER1.sam
OUTFILE=1.$ALIGNER1.accuracy

# writing filename to output
echo -e -n 'input_file\t' > $OUTFILE
echo $FILENAME >> $OUTFILE

# writing total reads
echo -e -n '\ntotal_reads\t' >> $OUTFILE
samtools view -c $FILENAME | head -c -1 >> $OUTFILE
echo -e -n '\ntotal_primary_reads\t' >> $OUTFILE
samtools view -cF256 $FILENAME | head -c -1 >> $OUTFILE

################################### collecting and writing read number diagnostics
echo -e -n '\n\nCM-S\t' >> $OUTFILE
correct_nuc_s=$(samtools view $FILENAME | awk -v x=$THRESHOLD '{if ($5 >= x){print $0}}' | awk '{if ($3 != "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
correct_mt_s=$(samtools view $FILENAME | awk -v x=$THRESHOLD '{if ($5 >= x){print $0}}' | awk '{if ($3 == "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
CM_S=$(( $correct_nuc_s + $correct_mt_s ))
echo -n $CM_S >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
echo -n $correct_nuc_s >> $OUTFILE
echo -e -n ')\t(' >> $OUTFILE
echo -n $correct_mt_s >> $OUTFILE
echo ')' >> $OUTFILE
###################################
echo -e -n 'CM-R\t' >> $OUTFILE
correct_nuc_r=$(samtools view -F256 $FILENAME | awk -v x=$THRESHOLD '{if ($5 >= x){print $0}}' | awk '{if ($3 != "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
correct_mt_r=$(samtools view -F256 $FILENAME | awk -v x=$THRESHOLD '{if ($5 >= x){print $0}}' | awk '{if ($3 == "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
CM_R=$(( $correct_nuc_r + $correct_mt_r ))
echo -n $CM_R >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
echo -n $correct_nuc_r >> $OUTFILE
echo -e -n ')\t(' >> $OUTFILE
echo -n $correct_mt_r >> $OUTFILE
echo ')' >> $OUTFILE
###################################
echo -e -n 'IM-S\t' >> $OUTFILE
incorrect_nuc_s=$(samtools view $FILENAME | awk -v x=$THRESHOLD '{if ($5 >= x){print $0}}' | awk '{if ($3 == "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
incorrect_mt_s=$(samtools view $FILENAME | awk -v x=$THRESHOLD '{if ($5 >= x){print $0}}' | awk '{if ($3 != "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
IM_S=$(( $incorrect_nuc_s + $incorrect_mt_s ))
echo -n $IM_S >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
echo -n $incorrect_nuc_s >> $OUTFILE
echo -e -n ')\t(' >> $OUTFILE
echo -n $incorrect_mt_s >> $OUTFILE
echo ')' >> $OUTFILE
###################################
echo -e -n 'IM-R\t' >> $OUTFILE
incorrect_nuc_r=$(samtools view -F256 $FILENAME | awk -v x=$THRESHOLD '{if ($5 >= x){print $0}}' | awk '{if ($3 == "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
incorrect_mt_r=$(samtools view -F256 $FILENAME | awk -v x=$THRESHOLD '{if ($5 >= x){print $0}}' | awk '{if ($3 != "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
IM_R=$(( $incorrect_nuc_r + $incorrect_mt_r ))
echo -n $IM_R >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
echo -n $incorrect_nuc_r >> $OUTFILE
echo -e -n ')\t(' >> $OUTFILE
echo -n $incorrect_mt_r >> $OUTFILE
echo ')' >> $OUTFILE
###################################
echo -e -n 'UM+Q\t' >> $OUTFILE
unmapped=$(samtools view -c -f4 $FILENAME)
qual=$(samtools view $FILENAME | awk -v x=$THRESHOLD '{if ($5 < x){print $0}}' | wc -l)
union=$(samtools view -f4 $FILENAME | awk -v x=$THRESHOLD '{if ($5 < x){print $0}}' | wc -l)
UM=$(( $unmapped+$qual-$union ))
echo -n $UM >> $OUTFILE
echo -e '\t(-)\t(-)' >> $OUTFILE
###################################
echo -e -n 'total\t' >> $OUTFILE
total=$(( $CM_S+$IM_S+$UM ))
echo -n $total >> $OUTFILE
echo -e -n '\t(nt)\t(mt)' >> $OUTFILE
###################################

# calculating and writing accuracy
echo -e -n '\n\nstrict_accuracy\t' >> $OUTFILE
printf '%3.4f\n' $(bc <<< "scale=4; $CM_S/($CM_S+$IM_S)") >> $OUTFILE
echo -e -n 'relaxed_accuracy\t' >> $OUTFILE
printf '%3.4f\n' $(bc <<< "scale=4; $CM_R/($CM_R+$IM_R)") >> $OUTFILE

# calculating and writing used read ratio
echo -e -n '\nused_read_ratio\t' >> $OUTFILE
printf '%3.4f\n' $(bc <<< "scale=4; ($CM_S+$IM_S)/($CM_S+$IM_S+$UM)") >> $OUTFILE

# detailed numbers
echo -e -n '\nunmapped_nuclear_reads\t' >> $OUTFILE
UNR=$(samtools view -f4 $FILENAME | cut -d '_' -f1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
echo -e $UNR >> $OUTFILE
echo -e -n 'unmapped_mitochondrial_reads\t' >> $OUTFILE
UMR=$(samtools view -f4 $FILENAME | cut -d '_' -f1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
echo -e $UMR >> $OUTFILE
echo -e -n 'mapped_nuclear_reads\t' >> $OUTFILE
MNR=$(samtools view -F4 $FILENAME | cut -d '_' -f1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
echo -e $MNR >> $OUTFILE
echo -e -n 'mapped_mitochondrial_reads\t' >> $OUTFILE
MMR=$(samtools view -F4 $FILENAME | cut -d '_' -f1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
echo -e $MMR >> $OUTFILE
echo -e -n 'total\t' >> $OUTFILE
echo -e $(( $UNR + $UMR + $MNR + $MMR )) >> $OUTFILE

# raw numbers (doesn't take MAPQ threshold into account)
echo -e -n '\nraw_CM\t' >> $OUTFILE
raw_correct_nuc=$(samtools view -F4 $FILENAME | awk '{if ($3 != "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
raw_correct_mt=$(samtools view -F4 $FILENAME | awk '{if ($3 == "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
RCM=$(( $raw_correct_nuc + $raw_correct_mt ))
echo -n $RCM >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
echo -n $raw_correct_nuc >> $OUTFILE
echo -e -n ')\t(' >> $OUTFILE
echo -n $raw_correct_mt >> $OUTFILE
echo ')' >> $OUTFILE
###################################
echo -e -n 'raw_IM\t' >> $OUTFILE
raw_incorrect_nuc=$(samtools view -F4 $FILENAME | awk '{if ($3 == "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
raw_incorrect_mt=$(samtools view -F4 $FILENAME | awk '{if ($3 != "MT"){print $1}}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
RIM=$(( $raw_incorrect_nuc + $raw_incorrect_mt ))
echo -n $RIM >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
echo -n $raw_incorrect_nuc >> $OUTFILE
echo -e -n ')\t(' >> $OUTFILE
echo -n $raw_incorrect_mt >> $OUTFILE
echo ')' >> $OUTFILE
###################################
echo -e -n 'raw_UM\t' >> $OUTFILE
echo -e -n $unmapped >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
echo -n $UNR >> $OUTFILE
echo -e -n ')\t(' >> $OUTFILE
echo -n $UMR >> $OUTFILE
echo ')' >> $OUTFILE
###################################
echo -e -n 'total\t' >> $OUTFILE
total=$(( $RCM+$RIM+$unmapped ))
echo -n $total >> $OUTFILE
echo -e '\t(nt)\t(mt)' >> $OUTFILE

# numts detailed breakdown
samtools view -h -F4 $FILENAME | awk '{if($3 != "MT"){print $0}}' > nDNA.mapped.sam
samtools view -b nDNA.mapped.sam | bedtools bamtobed -i - > nDNA.mapped.bed
###################################
echo -e -n '\nnuc reads mapped to nucleus\t' >> $OUTFILE
NRMTN=$(cat nDNA.mapped.bed | awk '{print $4}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
echo -n $NRMTN >> $OUTFILE
echo -e -n '\nmito reads mapped to nucleus\t' >> $OUTFILE
MRMTN=$(cat nDNA.mapped.bed | awk '{print $4}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
echo -n $MRMTN >> $OUTFILE
###################################
echo -e -n '\nnuclear reads mapped to nucleus and intersecting with numts\t' >> $OUTFILE
nuc_inter_numts=$(bedtools intersect -wa -a nDNA.mapped.bed -b $NUCLEARNUMTS | awk '{print $4}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
echo -n $nuc_inter_numts >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
printf '%3.2f' $(bc <<< "scale=2; ($nuc_inter_numts / $NRMTN) * 100") >> $OUTFILE
echo -e -n '%)' >> $OUTFILE
###################################
echo -e -n '\nnuclear reads mapped to nucleus & not intersecting with numts\t' >> $OUTFILE
nuc_notinter_numts=$(bedtools intersect -wa -a nDNA.mapped.bed -b $NUCLEARNUMTS -v | awk '{print $4}' | cut -d '_' -f 1 - | awk '{if ($0 != "MT"){print $0}}' | wc -l)
echo -n $nuc_notinter_numts >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
printf '%3.2f' $(bc <<< "scale=2; ($nuc_notinter_numts / $NRMTN) * 100") >> $OUTFILE
echo -e -n '%)' >> $OUTFILE
###################################
echo -e -n '\nmito reads mapped to nucleus and intersecting with numts\t' >> $OUTFILE
mito_inter_numts=$(bedtools intersect -wa -a nDNA.mapped.bed -b $NUCLEARNUMTS | awk '{print $4}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
echo -n $mito_inter_numts >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
printf '%3.2f' $(bc <<< "scale=2; ($mito_inter_numts / $MRMTN) * 100") >> $OUTFILE
echo -e -n '%)' >> $OUTFILE
###################################
echo -e -n '\nmito reads mapped to nucleus & not intersecting with numts\t' >> $OUTFILE
mito_notinter_numts=$(bedtools intersect -wa -a nDNA.mapped.bed -b $NUCLEARNUMTS -v | awk '{print $4}' | cut -d '_' -f 1 - | awk '{if ($0 == "MT"){print $0}}' | wc -l)
echo -n $mito_notinter_numts >> $OUTFILE
echo -e -n '\t(' >> $OUTFILE
printf '%3.2f' $(bc <<< "scale=2; ($mito_notinter_numts / $MRMTN) * 100") >> $OUTFILE
echo -e -n '%)' >> $OUTFILE

# removing temporary files
echo -e '\nremoving files.........'
rm nDNA.mapped.bed nDNA.mapped.sam

echo -e 'done\n'

done
