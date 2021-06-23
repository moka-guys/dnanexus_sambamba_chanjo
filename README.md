# dnanexus_sambamba_chanjo - v1.11

## What does this app do?
This app utilises Chanjo and Sambamba to calculates coverage.

### Sambamba
The `sambamba depth region` function is used with the following arguments:

* -L bedfile    -   Only counts regions stated in the bedfile
* -T n -   The required read depth (integer). For example, WES requires a minimum read depth of 20X, custom panels 30X whilst Oncology samples require a much higher read depth.
* -t `nproc`    -   Utilise multiple threads - uses the total number of threads available.
* -m    -   Does not count overlapping mate reads more than once.
* -F "mapping_quality >= 20"    - Uses the BAM flag mapping quality to only count bases with a mapping quality >=20
* -q    -   Min base quality for that base to be counted.

The sambamba output records the number of bases (that meet the parameters set) within each region of the bed file which have the required read depth.

This output is parsed by a python script to produce an output detailing any exons that are not covered 100% at the required read depth.

### Chanjo
The sambamba output is used by Chanjo to calculate the % of bases covered at the required read depth at a gene level.

The chanjo database is set up (`chanjo init -a; chanjo db setup`) and the sambamba file linked (`chanjo  link sambamba_output.bed`) and loaded (`chanjo load sambamba_output.bed`).

The bed file is used to extract a unique list of gene symbols and the coverage for each gene is then calculated (`chanjo calculate gene $gene >> chanjo_out.json`)

This output is then parsed by a custom python script to generate a file that can be downloaded, and loaded into MOKA to generate clinical coverage reports.


## What data are required for this app to run?
1. BAM file. This BAM file should be the same as used for variant calling, following all preprocessing.
2. The associated BAM file index (.bai)
3. A bed file in the sambamba format. This file is created by MokaBED and named in format Pan123sambamba.bed. 
4. The minimum read depth/coverage required (integer, default = 150). 
5. Count overlapping mate reads once? This True/False boolean denotes whether to apply the -m flag. True counts overlapping mate reads once (default) whereas False counts each overlapping mate read seperately.
6. Min base quality to be included (default = 25)
7. Min mapping quality score to be included (default = 20)
8. Exclude duplicate reads? This True/False boolean denotes whether to filter out reads that have failed quality control `- F not failed_quality_control`. True excludes reads that have failed quality control, false includes them.
9. There is an `additional_options` free form text box for passing additional flags to Sambamba, this excludes arguments passed as a text string to the -F filter flag.

## What does this app output?
This app produces five outputs:

1. The raw sambamba output file
2. The parsed sambamba output that identifies exons covered <100% at the stated coverage
3. The raw chanjo output
4. The parsed chanjo output that can be loaded into MOKA.
5. The chanjo.yaml file which provides the settings to chanjo

All files are output to a folder 'coverage'
the chanjo.yaml file is saved within a subfolder called 'chanjo_yaml' and the raw sambamba and chanjo files are output to a subfolder called 'raw_output'


## Created by
This app was created within the Viapath Genome Informatics section
