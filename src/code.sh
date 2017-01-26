#!/bin/bash
#

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#######################################################################
#
# Download inputs
#
#######################################################################

dx download "$sambamba_bed" -o bedfile
dx download "$bamfile"
dx download "$bam_index" 

# add miniconda to path to ensure correct version of python is used.
PATH=/home/dnanexus/miniconda2/bin:$PATH

#######################################################################
# run sambama
#######################################################################

# any errors with not being able to find the bedfile or reference file - check the bedfile for headers (and remove any!)

# if exome need 20x coverage
if [[ $bamfile_prefix == *"WES"* ]]; then
sambamba depth region -L bedfile -T 20 -m -F "mapping_quality >= 20" $bamfile_prefix.bam > sambamba_output.bed
printf "database: coverage.sqlite3\nsambamba:\n  cov_treshold:\n  - 20" > /home/dnanexus/chanjo.yaml
else 
#custom panells are at 30X
sambamba depth region -L bedfile -T 30 -m -F "mapping_quality >= 20" $bamfile_prefix.bam > sambamba_output.bed
printf "database: coverage.sqlite3\nsambamba:\n  cov_treshold:\n  - 30" > /home/dnanexus/chanjo.yaml
fi

head /home/dnanexus/chanjo.yaml
#######################################################################
# run chanjo
#######################################################################

#initiate database
chanjo init -a
chanjo db setup
# link bedfile
chanjo  link sambamba_output.bed
#load bedfile
chanjo  load sambamba_output.bed

# create a list of unique gene symbols
LIST=$(cat -n sambamba_output.bed | awk -F '\t' '$1>1 {print $9}' | sort | uniq)

#loop through list and get chanjo output
for i in $LIST
do
	chanjo calculate gene $i >> chanjo_out.json
done

#######################################################################
# clean up chanjo outputs using python script
#######################################################################
#python script parses the chanjo_out.json and sambamba_out, creating a chanjo_out.csv and a file
python read_chanjo.py

#######################################################################
# Upload results
#######################################################################

# make output folder
mkdir -p ~/out/chanjo_output/

# move and rename the chanjo out.json and sambamba_out
mv chanjo_out.json ~/out/chanjo_output/$bamfile_prefix.chanjo_out_json
mv sambamba_output.bed ~/out/chanjo_output/$bamfile_prefix.sambamba_output.bed
mv /home/dnanexus/chanjo.yaml ~/out/chanjo_output/chanjo.yaml

#move and rename the output of the python script
mv chanjo_out.txt ~/out/chanjo_output/$bamfile_prefix.chanjo_txt
mv exon_level.txt ~/out/chanjo_output/$bamfile_prefix.exon_level.txt

dx-upload-all-outputs --parallel
