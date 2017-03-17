# dnanexus_sambamba_chanjo

## What does this app do?
This app utilises Chanjo and Sambamba to calculates coverage.

The script reads the name of the BAM file to determine if it is a WES sample or a custom panel, and therefore if coverage should be calculated at 20X or 30X.
####Sambamba
The sambamba depth region function is used. 
We provide the following arguments:

* -L bedfile    -     Only counts regions stated in the bedfile

* -T 20         -      The required read depth of 20X (Custom panels has -T 30)

* -m            - Does not count overlapping mate reads more than once.

* -F "mapping_quality >= 20"  - Uses the BAM flag mapping quality to only count bases with a mapping quality >=20


The sambamba output records the number of bases (that meet the parameters set) within each region of the bed file which have the required read depth.

This output is parsed by a python script into a format that can be easily used to identify exons which are not covered 100% at >=30X and therefore require filling via sanger sequencing.

####Chanjo
The sambamba output is used by Chanjo to calculate the % of bases covered at the required read depth at a gene level. 

This output is then parsed by a custom python script to generate a file that can be downloaded, and loaded into MOKA to generate clinical coverage reports.


## What data are required for this app to run?

This app requires three input files:

1.  A BAM file. Where possible this should be the refined bam

2. The BAM file index (.bai)

3. A bed file in the sambamba format. This file is created by MokaBED and named Pan100sambamba.bed. 



## What does this app output?
This app produces five outputs:

1. The raw sambamba output file

2. The parsed sambamba output that identifies exons covered <100% at 30X

3. The raw chanjo output

4. The parsed chanjo output that can be loaded into MOKA.

5. The chanjo.yaml file which provides the settings to chanjo


All files are output to a folder 'coverage'
the chanjo.yaml file is saved within a subfolder called 'chanjo_yaml' and the raw sambamba and chanjo files are output to a subfolder called 'raw_output'


## Created by
This app was created within the Viapath Genome Informatics section
