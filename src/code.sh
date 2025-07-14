#!/bin/bash

# The following line causes bash to exit at any point if there is any error
# and to output each line as it is executed -- useful for debugging
set -e -x -o pipefail

#######################################################################
# Install Python dependencies
#######################################################################
pip3 install reportlab

#######################################################################
# Download inputs
#######################################################################

dx download "$sambamba_bed" -o bedfile

# Get the actual BAM filename and download with original name
bamfile_name=$(dx describe "$bamfile" --json | jq -r '.name')
dx download "$bamfile" -o "$bamfile_name"
dx download "$bam_index" -o "${bamfile_name}.bai"

#######################################################################
# Extract BAM filename prefix for outputs
#######################################################################

bamfile_prefix=$(basename "$bamfile_name" .bam)

#######################################################################
# Run coverage analysis using the Python script
#######################################################################

# The Python script handles sambamba execution and report generation
# Pass all the DNAnexus inputs as arguments to the Python script

# Build the command with conditional arguments
coverage_cmd="python3 /home/dnanexus/coverage.py --bam-file \"${bamfile_name}\" --bed-file bedfile --coverage-level ${coverage_level} --min-mapping-qual ${min_mapping_qual} --min-base-qual ${min_base_qual} --additional-filter \"${additional_filter_commands}\" --threads $(nproc) --output-dir coverage_results --panel-filter --pdf-report"

# Add conditional flags based on boolean inputs
if [[ "$exclude_failed_quality_control" == true ]]; then
    coverage_cmd="$coverage_cmd --exclude-failed-qc"
fi

if [[ "$exclude_duplicate_reads" == true ]]; then
    coverage_cmd="$coverage_cmd --exclude-duplicates"
fi

if [[ "$merge_overlapping_mate_reads" == true ]]; then
    coverage_cmd="$coverage_cmd --merge-overlapping-mates"
fi

# Add additional sambamba flags if provided
if [ -n "$additional_sambamba_flags" ]; then
    coverage_cmd="$coverage_cmd --additional-sambamba-flags \"${additional_sambamba_flags}\""
fi

# Execute the coverage analysis
eval $coverage_cmd

#######################################################################
# Upload results
#######################################################################

# Create output directories
mkdir -p ~/out/sambamba_output/coverage/
mkdir -p ~/out/gene_level_report/coverage/
mkdir -p ~/out/exon_level_report/coverage/
mkdir -p ~/out/pdf_report/coverage/

# Move and rename the outputs from the Python script
mv coverage_results/sambamba_output.bed ~/out/sambamba_output/coverage/${bamfile_prefix}.sambamba_output.bed
mv coverage_results/gene_level.txt ~/out/gene_level_report/coverage/${bamfile_prefix}.gene_level.txt
mv coverage_results/exon_level.txt ~/out/exon_level_report/coverage/${bamfile_prefix}.exon_level.txt

# Move PDF report if it exists
if [ -f "coverage_results/${bamfile_prefix}_coverage_report.pdf" ]; then
    mv coverage_results/${bamfile_prefix}_coverage_report.pdf ~/out/pdf_report/coverage/${bamfile_prefix}_coverage_report.pdf
else
    echo "PDF report not generated, skipping PDF output"
fi

# Upload all outputs
dx-upload-all-outputs --parallel