#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#######################################################################
# Download inputs
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
# Use the coverage_level input to specify the coverage to be reported
# Use sambaba depth to append coverage to sambamba_output.bed 
# -L is the bedfile to define the regions that coverage must be calculated for
# -T is the minimum coverage required for the amplicon
# -t is the number of threads available
# -m does not count overlapping mate reads more than once
# -F allows filtering using other BAM info eg mapping quality
# check if -m flag is required

opts="-T $coverage_level -q $min_base_qual"

if [[ "$merge_overlapping_mate_reads" == true ]]; then
	opts="$opts -m"
fi

sambamba depth region -L bedfile -t "$(nproc)"  "$opts" -F "$min_mapping_qual >= 20 and not duplicate and not failed_quality_control" "$bamfile_prefix".bam > sambamba_output.bed




# chanjo can be modified using a yaml config file. The threshold level can be specified here using coverage_level input.
printf "database: coverage.sqlite3\nsambamba:\n  cov_treshold:\n  - $coverage_level" > /home/dnanexus/chanjo.yaml


#head sambamba_output.bed
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
for gene in $LIST
do
	chanjo calculate gene $gene >> chanjo_out.json
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
mkdir -p ~/out/chanjo_raw_output/coverage/raw_output/
mkdir -p ~/out/chanjo_yaml/coverage/chanjo_yaml/
mkdir -p ~/out/chanjo_output_to_report/coverage/

# move and rename the chanjo out.json and sambamba_out
mv chanjo_out.json ~/out/chanjo_raw_output/coverage/raw_output/$bamfile_prefix.chanjo_out_json
mv sambamba_output.bed ~/out/chanjo_raw_output/coverage/raw_output/$bamfile_prefix.sambamba_output.bed
mv /home/dnanexus/chanjo.yaml ~/out/chanjo_yaml/coverage/chanjo_yaml/$bamfile_prefix.chanjo.yaml

#move and rename the output of the python script
mv chanjo_out.txt ~/out/chanjo_output_to_report/coverage/$bamfile_prefix.chanjo_txt
mv exon_level.txt ~/out/chanjo_output_to_report/coverage/$bamfile_prefix.exon_level.txt


dx-upload-all-outputs --parallel