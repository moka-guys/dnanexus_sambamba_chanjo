#!/bin/bash
#

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#
# Fetch bed files and bam files
#
mark-section "Downloading inputs"
dx download "$sambamba_bed" -o bedfile
dx download "$bamfile"
dx download "$bam_index" 

#
# add miniconda to path to ensure correct version of python is used.
#
PATH=/home/dnanexus/miniconda2/bin:$PATH

#
# run sambama
#
mark-section "running sambamba"
sambamba depth region -L bedfile -T 20 -F "mapping_quality >= 20" $bamfile_prefix.bam > $bamfile_prefix.sambamba_output.bed

#
# run chanjo
#
mark-section "running chanjo"
#initiate database
chanjo init -a 
chanjo db setup
# link bedfile
chanjo  link $bamfile_prefix.sambamba_output.bed
#load bedfile
chanjo  load $bamfile_prefix.sambamba_output.bed

# create a list of unique gene symbols
LIST=$(cat -n $bamfile_prefix.sambamba_output.bed | awk -F '\t' '$1>1 {print $9}' | sort | uniq)

#loop through list and get chanjo output
for i in $LIST
do
	chanjo calculate gene $i >> chanjo_out.json
done

mark-section "running python script to parse chanjo output"
#python script parses the chanjo_out.json and creates a chanjo_out.csv file
python read_chanjo.py

#
# Upload result
#
mark-section "uploading results"
# make output folder
mkdir -p ~/out/chanjo_output/
#move and rename the output of the python script
mv chanjo_out.csv ~/out/chanjo_output/$bamfile_prefix.chanjo_out
dx-upload-all-outputs --parallel
