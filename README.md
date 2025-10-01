# dnanexus_sambamba_chanjo - v2.0

## What does this app do?

This applet calculates sequencing coverage statistics using Sambamba and a custom Python script.

### Sambamba

The process is a simple, two-step pipeline:

- Sambamba: The sambamba depth region command is run first. It uses a BED file to calculate coverage statistics (mean coverage, percentage of bases above a threshold, etc.) for each genomic region.

- Python Parser: A custom Python script (coverage.py) then reads the raw sambamba output file to generate two final, user-friendly reports: a gene-level summary and an exon-level "report-by-exception".


The sambamba depth region function is used with the following arguments:

    -L: The BED file defining the regions of interest.

    -T: The required read depth threshold (e.g., 30 for 30X).

    -t: The number of threads to use for processing.

    -m: An optional flag to count overlapping mate pairs only once.

    -F: A powerful filter string to include or exclude reads based on BAM flags (e.g., mapping quality, duplicates).

    --min-base-quality: The minimum base quality score for a base to be included in coverage calculations.


## What data are required for this app to run?

1. BAM file. This BAM file should be the same as used for variant calling, following all preprocessing.
2. The associated BAM file index (.bai)
3. A bed file in the sambamba format. This file is created by MokaBED and named in format Pan123sambamba.bed.
4. The minimum read depth/coverage required (integer).
5. Count overlapping mate reads once? This True/False boolean denotes whether to apply the -m flag. True counts overlapping mate reads once (default) whereas False counts each overlapping mate read seperately.
6. Min base quality to be included (default = 25)
7. Min mapping quality score to be included (default = 20)
8. Exclude duplicate reads? This True/False boolean denotes whether to filter out reads that have failed quality control `- F not failed_quality_control`. True excludes reads that have failed quality control, false includes them.
9. There is an `additional_sambamba_flags` free form text box for passing additional flags to Sambamba, this excludes arguments passed as a text string to the -F filter flag which are passed using `additional_filter_commands`.
10. There is an `additional_filter_commands` free form text box for passing additional arguments as a text string to the Sambamba -F filter flag.  Commands are appended to the existing filter.

## What does this app output?

All outputs are saved into a folder 'coverage'. Files that are not needed routinely saved into subfolders.
This app produces three outputs:

1. [sample_name].gene_level.txt: A high-level summary of coverage for every gene. This report includes:

    gene_symbol: The name of the gene.

    mean_coverage: The mean coverage across all exons in the gene, weighted by exon length.

    average_completeness: The percentage of bases that met the coverage threshold, also weighted by exon length.

    exon_count: The total number of exons analyzed for the gene.

    total_length_bp: The total length of the analyzed exons in base pairs.

2. [sample_name].exon_level.txt: A "report-by-exception" that makes it easy to find regions with poor coverage. This file only lists exons that are covered less than 100% at the specified depth. It is separated into "CODING" and "UTR" sections for clarity.

3. [sample_name].sambamba_output.bed: The raw, complete output from the sambamba depth region command. This file is useful for manual review or debugging.

## Running via the command line

The `dx-toolkit` can be used to run the app from the command line using the format below:

```bash
dx run 001_ToolsReferenceData:/Apps/chanjo_sambamba_coverage_v2.0 \
--name "congenicaV2.0_test_standard settings" \
-i sambamba_bed=project-PXH0qBZHXlfT2JYf7aMcHMxQk:file-cgEB1sJsCtyIzNnD6rju57Nhe \
-i bamfile=project-PXH0qBZHXlfT2JYf7aMcHMxQk:file-Bk5bbrm6YuQGjAIupujWk0FjO \
-i bam_index=project-PXH0qBZHXlfT2JYf7aMcHMxQk:file-BppUJ6p6cJGPS2WWEIWqzMBOc \
-i merge_overlapping_mate_reads=true \
-i exclude_failed_quality_control=true \
-i exclude_duplicate_reads=true \
-i coverage_level="150" \
-i additional_sambamba_flags="--annotate" \
-i additional_filter_commands="not mate_is_unmapped and first_of_pair" \
-i min_base_qual="25" \
-i min_mapping_qual="20" \
--brief
```

## Testing

There are two modules utilising pytest which can be invoked by running the following commands from within the app's folder:

```bash
# Run all modules
pytest --setup-show
# Individual modules can be run using the format:
pytest /src/test_inputs.py --setup-show
```

* `test_inputs.py` - this tests that various combinations of input parameters don't break the app.
* `test_outputs.py` - this tests data run through a previous version of the app for MokaWES, MokaPipe, and MokaAmp and checks that the output is similar.
### NOTE: Currently test_outputs.py tests that Sambamba will run with the output of multiple pipelines, but does not automatically compare that output.  This can be done manually by comparing the produced files to those from the original.
## Created by

This app was created within the Synnovis Bioinformatics team
