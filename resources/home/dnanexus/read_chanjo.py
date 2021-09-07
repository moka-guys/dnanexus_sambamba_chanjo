#!/usr/bin/python

'''
Created on 14 Dec 2016

@author: aled
'''
import sys
import subprocess 
import os

class read_chanjo():
    def __init__(self):
        self.path="/home/dnanexus/"
        #self.path="/home/aled/Documents/210906_gene_level_coverage/exome/"
        self.json_file = self.path+"chanjo_out.json"
        self.output_file=self.path+"chanjo_out.txt"
        self.sambamba_in=self.path+"sambamba_output.bed"
        self.exonlevel=self.path+"exon_level.txt"
        self.genelevel=self.path+"gene_level.txt"
        # Read the bash variable which contains the minimum coverage level (specified in the app's $coverage_level input)
        self.coverage_level = str(os.environ['coverage_level'].rstrip())


    def read_json(self):
        # open output file
        output=open(self.output_file, 'w')
        # loop through json
        with open(self.json_file,'r') as json:
            for line in json:
                #split line
                line_split= line.split(",")
                # capture gene, percent_bases_covered and average coverage
                gene=line_split[0].split(":")[1].replace(" {","").replace("\"","").replace("{","")
                percent_bases_covered=line_split[0].split(":")[3]
                average=line_split[1].split(":")[1].replace("}}","")
                #write to output
                output.write(gene+"\t"+percent_bases_covered+"\t"+average+"\n")    
        # close output file
        output.close()

    def exon_level(self):
        # open exon level coverage output
        output=open(self.exonlevel, 'w')
        # write header
        output.write("gene\ttranscript\tentrezID\tChr\tstart\tstop\taverage_coverage\tpercent_bases_covered\n")    
        
        #open sambamba output
        with open(self.sambamba_in,'r') as samb:
            # ignore header
            for line in samb:
                if not line.startswith("#"):
                    # split and capture gene name, entrezid, coordinates, average coverage and the percent_bases_covered
                    line_split= line.split("\t")
                    gene=line_split[6]
                    entrezid=line_split[7].rstrip()
                    coords=line_split[3]
                    meanCoverage=float(line_split[9])
                    percent_bases_covered=float(line_split[10])            
                    # if not covered 100% @ required depth 
                    if percent_bases_covered < 100.00:
                        # some bedfiles now have rsIDs in the gene column so are not in the format gene;transcript.
                        # if it's in format gene;transcript
                        if len(gene.split(";")) == 2:
                            coverage_report_gene = gene.split(";")[0]
                            coverage_report_transcript = gene.split(";")[1]
                        # if it's an rsID
                        elif len(gene.split(";")) == 1 and gene.startswith("rs"):
                            coverage_report_gene = gene
                            coverage_report_transcript = gene
                        else:
                            raise AssertionError("unable to find the gene symbol")
                        # write to file
                        output.write(coverage_report_gene+"\t"+coverage_report_transcript+"\t"+entrezid+"\t"+coords.split("-")[0]+"\t"+coords.split("-")[1]+"\t"+coords.split("-")[2]+"\t"+str(meanCoverage)+"\t"+str(percent_bases_covered)+"\n")    
            else:
                # write line to end report
                output.write("Any exons not mentioned above are covered 100% at " + self.coverage_level + "X")
        
        # close output file
        output.close()        

    def gene_level(self):
        """
        A human readable gene level report is required to describe the % of each gene covered at the given read depth
        This takes the chanjo json file and uses the entrezgeneid to search the sambamba BED file to pull out the gene symbol
        """
        # open exon level coverage output
        output=open(self.genelevel, 'w')
        # write header
        output.write("\t".join(["gene symbol","percent_bases_covered at %s\n" % (self.coverage_level + "X")]))
        # loop through newly created chanjo_out.txt - this has 3 columns- capture gene, percent_bases_covered and average coverage
        with open(self.output_file,'r') as json:
            for line in json:
                # we want to capture the gene_symbol - set as empty variable for each gene/line in the input file
                gene_symbol = ""
                entrezgeneid,percent_bases_covered,average = line.split("\t")
                # open and loop through sambamba.bed to match on entrezgeneid and then pull out corresponding gene symbol
                with open(self.sambamba_in) as sambamba_bed:
                    for line in sambamba_bed.readlines():
                        chrom,chromStart,chromEnd,F3,F4,F5,F6,F7,readCount,meanCoverage,percentage600,sampleName = line.split("\t")
                        # the column containing gene symbol is in format genesymbol;transcript
                        # however the sambamba/bed includes some other regions, such as SNV control sites 
                        # these have the entrezgeneid set to 0 and the column F6 in format ";rs123"
                        # There are also clinically relevant SNVs whitelisted so we don't see incidental findings elsewhere in the gene eg CHEK2. 
                        # These have the correct entrezgeneid present but unfortunately the gene symbol is not in column F6 (";rs123")
                        # Therefore we need to hard code this gene symbol (until the bed file is corrected)
                        if str(entrezgeneid) == "11200":
                            gene_symbol = "CHEK2"
                        # for all other entrezgeneids look to match columns but ignore control sites where entrezgeneid = 0.
                        elif str(entrezgeneid) == str(F7) and str(entrezgeneid) != "0":
                            gene_symbol = F6.split(";")[0]
                # if the genesymbol has been identified
                if gene_symbol:
                    # the percent bases covered column appears to have a leading space. not wanted to change this incase it breaks downstream processes. use lstrip to remove
                    output.write("\t".join([gene_symbol,percent_bases_covered.lstrip()+"\n"]))
                # # used for testing
                else:
                    print "unable to find a gene symbol for entrezgeneid %s" % (entrezgeneid)
        # close output file
        output.close()
if __name__ == '__main__':
    a=read_chanjo()
    a.read_json()
    a.exon_level()
    a.gene_level()