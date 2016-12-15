#!/bin/bash
#

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch reads
#
dx download "$sambamba_bed" -o bedfile
dx download "$bamfile"
dx download "$bam_index" 

#
# set python path
#
PATH=/home/dnanexus/miniconda2/bin:$PATH
ls
#
# run sambama
#
#bedtools sort -i bedfile > $sambamba_bed_prefix.sorted.bed
sambamba depth region -L bedfile -T 20 -F "mapping_quality >= 20" $bamfile_prefix.bam > $bamfile_prefix.sambamba_output.bed

#
# run chanjo
#
chanjo init -a 
chanjo db setup
chanjo  link $bamfile_prefix.sambamba_output.bed
chanjo  load $bamfile_prefix.sambamba_output.bed
LIST=$(cat -n $bamfile_prefix.sambamba_output.bed | awk -F '\t' '$1>1 {print $9}' | sort | uniq)
for i in $LIST
do
	chanjo calculate gene $i >> chanjo_out.json
done
pwd
ls
python read_chanjo.py

#
# Upload result
#
mkdir -p ~/out/chanjo_output/
mv chanjo_out.csv ~/out/chanjo_output/$bamfile_prefix.chanjo_out.csv
dx-upload-all-outputs --parallel
