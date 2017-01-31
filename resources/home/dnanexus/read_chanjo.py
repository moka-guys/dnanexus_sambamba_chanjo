#!/usr/bin/python

'''
Created on 14 Dec 2016

@author: aled
'''
import sys

class read_chanjo():
    def __init__(self):
        self.path="/home/dnanexus/"
        self.json_file = self.path+"chanjo_out.json"               
        self.output_file=self.path+"chanjo_out.txt"
        self.sambamba_in=self.path+"sambamba_output.bed"
        self.exonlevel=self.path+"exon_level.txt"
    
        #self.json_file = "/home/aled/Documents/DNA_Nexus_app_github/chanjo2/4221_S5_L001_R1_001.chanjo_out.json" 
        #self.output_file="/home/aled/Documents/DNA_Nexus_app_github/chanjo2/4221_S5_L001_R1_001.chanjo_out.csv"

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
        #pen sambamba output
        with open(self.sambamba_in,'r') as samb:
            # ignore header
            for line in samb:
                if line.startswith("#"):
                    pass
                else:
                    # split and capture gene name, coordinates and the percent_bases_covered
                    line_split= line.split("\t")
                    gene=line_split[6]
                    coords=line_split[3]
                    percent_bases_covered=line_split[10]
                    # write to file
                    output.write(gene+"\t"+coords+"\t"+percent_bases_covered+"\n")    
        
        # close output file
        output.close()        


if __name__ == '__main__':
    a=read_chanjo()
    a.read_json()
    a.exon_level()
    