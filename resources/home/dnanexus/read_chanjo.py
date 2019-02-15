#!/usr/bin/python

'''
Created on 14 Dec 2016

@author: aled
'''
import sys
import subprocess 

class read_chanjo():
    def __init__(self):
        self.path="/home/dnanexus/"
        self.json_file = self.path+"chanjo_out.json"               
        self.output_file=self.path+"chanjo_out.txt"
        self.sambamba_in=self.path+"sambamba_output.bed"
        self.exonlevel=self.path+"exon_level.txt"
    
        #self.json_file = "/home/aled/Documents/DNA_Nexus_app_github/chanjo2/4221_S5_L001_R1_001.chanjo_out.json" 
        #self.output_file="/home/aled/Documents/DNA_Nexus_app_github/chanjo2/4221_S5_L001_R1_001.chanjo_out.csv"
        self.coverage_level = "X"

    def get_coverage_level(self):
        cmd = "echo $coverage_level"
        proc = subprocess.Popen([cmd], stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True)
        # Capture the streams
        (out, err) = proc.communicate()
        print "out"+out
        print "err"+err
        self.coverage_level = out.rstrip() + self.coverage_level
        print self.coverage_level

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
        
        # only report regions which are not covered at 100%
        count=0

        #open sambamba output
        with open(self.sambamba_in,'r') as samb:
            # ignore header
            for line in samb:
                if line.startswith("#"):
                    #pass
                    output.write("gene\ttranscript\tentrezID\tChr\tstart\tstop\taverage_coverage\tpercent_bases_covered\n")    
                else:
                    # split and capture gene name, coordinates and the percent_bases_covered
                    line_split= line.split("\t")
                    gene=line_split[6]
                    entrezid=line_split[7].rstrip()
                    coords=line_split[3]
                    meanCoverage=float(line_split[9])
                    percent_bases_covered=float(line_split[10])            
                    if percent_bases_covered < 100.00:
                        # write to file
                        output.write(gene.split(";")[0]+"\t"+gene.split(";")[1]+"\t"+entrezid+"\t"+coords.split("-")[0]+"\t"+coords.split("-")[1]+"\t"+coords.split("-")[2]+"\t"+str(meanCoverage)+"\t"+str(percent_bases_covered)+"\n")    
                        count += 1
            else:
                output.write("Any exons not mentioned above are covered 100% at " + self.coverage_level)
        
        # close output file
        output.close()        


if __name__ == '__main__':
    a=read_chanjo()
    a.get_coverage_level()
    a.read_json()
    a.exon_level()
    