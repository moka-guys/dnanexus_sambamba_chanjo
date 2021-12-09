<<<<<<< HEAD
# dnanexus_sambamba_chanjo - v1.13
=======
# dnanexus_sambamba_chanjo - v1.12
>>>>>>> a57c79dfddfb157b5f78c08162f72a64740d61bb

## What does this app do?

This app utilises Chanjo and Sambamba to calculates coverage.

### Sambamba

The `sambamba depth region` function is used with the following arguments:

* -L bedfile    -   Only counts regions stated in the bedfile
* -T n -   The required read depth (integer). For example, WES requires a minimum read depth of 20X, custom panels 30X whilst Oncology samples require a much higher read depth.
* -t $(nproc)    -   Utilise multiple threads - uses the total number of threads available.
* -m    -   Does not count overlapping mate reads more than once.
* -F "mapping_quality 20 and not failed_quality_control and not duplicate"    - DEFAULT: Uses the BAM flag mapping quality to only count bases with a mapping quality >=20, which have passed QC, and are not duplicates. Syntax can be found [here](https://github.com/biod/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax)
* -q    -   Min base quality for that base to be counted.

The sambamba output records the number of bases (that meet the parameters set) within each region of the bed file which have the required read depth.

This output is parsed by a python script to produce an output detailing any exons that are not covered 100% at the required read depth.

### Chanjo

The sambamba output is used by Chanjo to calculate the % of bases covered at the required read depth at a gene level.

The chanjo database is set up (`chanjo init -a; chanjo db setup`) and the sambamba file linked (`chanjo  link sambamba_output.bed`) and loaded (`chanjo load sambamba_output.bed`).

The bed file is used to extract a unique list of gene symbols and the coverage for each gene is then calculated (`chanjo calculate gene $gene >> chanjo_out.json`)
### Read_chanjo.py
A python script parses the sambamba and chanjo outputs and generates three files in formats that can be used downstream (see below)
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
This app produces six outputs.

1. Exon level coverage (*exon_level.txt). This consists of any exons which are covered <100% at the stated coverage (coverage/)
2. gene wide coverage for moka (*chanjo_txt). This takes the chanjo output and reformats it so coverage can be inserted into Moka (WES samples). For each gene the entrezgene id, % covered at the the given X and average coverage in a tab seperated list (coverage/)
3. Gene level coverage report (*gene_level.txt). Introduced in v1.12. This is a further modification of the chanjo output, but including the human readable gene symbol for the analyst (coverage/)
4. The raw chanjo output (*chanjo_out.json). The chanjo output in JSON format (coverage/raw_output)
5. The raw sambamba output file(*sambamba_output.bed). The sambamba BED file that is input to the app with columns describing coverage added to the end of each row (coverage/raw_output)
6. The chanjo.yaml file (*chanjo.yaml). YAML file which provides the settings to chanjo (coverage/chanjo_yaml)

## Running via the command line

The `dx-toolkit` can be used to run the app from the command line using the format below:

```bash
dx run 001_ToolsReferenceData:/Apps/chanjo_sambamba_coverage_v1.13 \
--name "congenicaV1.13_test_standard settings" \
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

There are two modules utilizing pytest which can be invoked by running the following commands from within the app's folder:

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

This app was created within the Viapath Genome Informatics section
