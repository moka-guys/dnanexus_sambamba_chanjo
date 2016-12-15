#!/usr/bin/python

'''
Created on 14 Dec 2016

@author: aled
'''

class read_chanjo():
    def __init__(self):
        self.json_file = "/home/dnanexus/chanjo_out.json"               
        self.output_file="/home/dnanexus/chanjo_out.csv"
        #self.json_file = "/home/aled/Documents/DNA_Nexus_app_github/chanjo2/4221_S5_L001_R1_001.chanjo_out.json" 
        #self.output_file="/home/aled/Documents/DNA_Nexus_app_github/chanjo2/4221_S5_L001_R1_001.chanjo_out.csv"

    def read_json(self):
        output=open(self.output_file, 'w')

        with open(self.json_file,'r') as json:
            for line in json:
                line_split= line.split(",")
                gene=line_split[0].split(":")[1].replace(" {","").replace("\"","").replace("{","")
                twentyx=line_split[0].split(":")[3]
                average=line_split[1].split(":")[1].replace("}}","")
                output.write(gene+"\t"+twentyx+"\t"+average+"\n")    
        output.close()


if __name__ == '__main__':
    a=read_chanjo()
    a.read_json()
    