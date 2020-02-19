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
        self.json_file = self.path+"chanjo_out.json"               
        self.output_file=self.path+"chanjo_out.txt"
        self.sambamba_in=self.path+"sambamba_output.bed"
        self.exonlevel=self.path+"exon_level.txt"
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


if __name__ == '__main__':
    a=read_chanjo()
    a.read_json()
    a.exon_level()
    