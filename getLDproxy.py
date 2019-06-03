#!/usr/bin/python3

#Script by Brooke Wolford
#June 3, 2019

import gzip
import sys
import re
import argparse
import pysam
from collections import namedtuple
from contextlib import closing
#from requests import HTTPSession
import requests

argparser = argparse.ArgumentParser('Finds proxy-SNPs given coodinates or rsID')
argparser.add_argument("--in_rs",type=str,help="Inputt file with new line separated rsIDs")
argparser.add_argument("--in_coord",type=str,help="Input file with new line separated coordinates X:XXX")
argparser.add_argument("--in_vcf",required=True,type=str,help="VCF file from dbSNP with coordinates and rsIDs. .tbi file required in the same directory.")
argparser.add_argument("--proxy",action='store_true',help="Look up proxy information for all SNPs provided")

def open_zip(f):
    if ".gz" in f:
        command=gzip.open(f,"rt")
        print >> sys.stderrr, "Opening gzipped file %s\n" % f
    elif f == "-":
        sys.exit("Cannot read file from stdin\n")
    else:
        command=open(f,"rt")
        print >> sys.stderr, "Opening file %s\n" % f
        return command
    
    def get_rsID(coord_file,vcf):
        
        return 0
    
    def get_coord(rs_file,vcf):
        command=open_zip(rs_file)
        with command as f:
            for line in f:
                ls=line.rstrip()
        query_string=''.join("'-i'%ID=@",rs_file,"'")        
        subprocess.call(["bcftools","query",query_string, "-f'%CHROM:%POS'",vcf])
        # subprocess.call([bcftools,"query",vcf,"-R",tmp.name,"-f","%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n","-o",outName])
        
        return 0
            
    def get_proxy():
        r=requests.get("https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=rs3&pop=MXL&r2_d=r2&token=5281fa13d3e3")
        print(r.text)
    #use Brooke's API token for  https://ldlink.nci.nih.gov/?tab=apiaccess
    #curl -k -X GET 'https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=rs3&pop=MXL&r2_d=r2&token=5281fa13d3e3'
#does it need coordinattes or rsids
                    
    return 0
            
if __name__ == '__main__':
    args = argparser.parse_args()
    print(args)
    
    if args.in_rs is not None and args.in_coord is not None:
        raise Exception('Supply rsID or coordinates but not both\n')
    elif args.in_rs is not None:
        get_coord(args.in_rs,args.in_vcf)
    elif args.in_coord is not None:
        get_rsID(args.in_rs,args.in_vcf)
    else:
        raise Exception('Supply either rsID or coordinate file\n')
