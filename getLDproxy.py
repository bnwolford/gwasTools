#!/usr/bin/python3

#Script by Brooke Wolford
#June 3, 2019

import gzip
import sys
import re
import argparse
import subprocess
from collections import namedtuple
from contextlib import closing
#from requests import HTTPSession
import requests
from tempfile import NamedTemporaryFile

### Written for Python 2.7

argparser = argparse.ArgumentParser('Finds proxy-SNPs given coodinates or rsID and checks their availability in a study-speciifc VCF file. Outputs *bcftools.txt which is the bcftools query output, *proxy.txt which is the results from the LDproxy query and *study.txt which is the snp id for snps of interest or their proxies that is in the study specific VCF')
argparser.add_argument("--in_rs",type=str,help="Inputt file with new line separated rsIDs")
argparser.add_argument("--in_coord",type=str,help="Input file with new line separated coordinates X:XXX")
argparser.add_argument("--dbsnp_vcf",required=True,type=str,help="VCF file from dbSNP with coordinates and rsIDs. .tbi file required in the same directory.")
argparser.add_argument("--population",help="Population code for r2 info [default=CEU]",default="CEU",type=str)
argparser.add_argument("--minrsq",help="Minimum rsq for printing variant",default=0.8, type=float)
argparser.add_argument("--token",type=str, default="5281fa13d3e3",help="Token, sign up here https://ldlink.nci.nih.gov/?tab=apiaccess",required=True)
argparser.add_argument("--convert",action='store_true',help="Convert rsid to coordinate or vice versa")
argparser.add_argument("--study_vcf",type=str,help="VCF with sites that are analyzed in HUNT",required=True)
argparser.add_argument("--output",type=str,help="Output prefix [default=LDproxy]",default="LDproxy")
argparser.add_argument("--bcftools",type=str,help="Path to bcftools [default='bcftools']",default="bcftools")

def open_zip(f):
    """open files"""
    if ".gz" in f:
        command=gzip.open(f,"rt")
        print >> sys.stderrr, "Opening gzipped file %s\n" % f
    elif f == "-":
        sys.exit("Cannot read file from stdin\n")
    else:
        command=open(f,"rt")
        print >> sys.stderr, "Opening file %s\n" % f
        return command
    
def read_coord(coord_file,vcf,convert):
    print >> sys.stderr, "This function not available yet"
    return 0
    
def read_rsid(rs_file,vcf,convert,token,pop,minrsq,out):
    """ Read all the rsIDs and look for LD buddies. Create a flat file of SNPs of interest and LD buddies if found. """
    if convert==True:
        #TO DO: convert rsID to coordinate
        print >> sys.stderr, "This function not available yet"
        return 0
    else:
        rsid_list=[] #initialize rsid list
        command=open_zip(rs_file)
        with command as f:
            for line in f:
                ls=line.rstrip()
                rsid_list.append(ls)
        count=0
        error_list=[] #initialize list of SNPs of interest with failed queries
        fn=".".join([out,"proxy.txt"]) #filename
        with open(fn, 'w') as output: #open file to write to
            
            for r in rsid_list: #for every SNP of interest
                query_list=get_proxy(token,pop,minrsq,r) #perform query for LD proxies
                if query_list!=20: #if no custom error code 20 from get_proxy
                    if count==0: #print header once for all SNPs of interest
                        output.write("\t".join(["\t".join(query_list[0]),"SNP_of_interest"]))
                        output.write("\n")
                        count=count+1
                    for proxy in query_list[1:]: #print 
                        output.write("\t".join(["\t".join(proxy),r]))
                        output.write("\n")
                else: #if error, save the snps of interest to output at the end 
                    error_list.append(r)
        output.close()
        print(len(error_list))
        
        #get coordinates for SNPs with rsID and no proxy so we can still search study for them
        noLD_out=".".join([out,"noLD.txt"]) #write a file with one rsID per line of the SNPs without LD buddies
        with open(noLD_out,"w") as noLD:
            for err in error_list:
                noLD.write(err)
                noLD.write("\n")
        noLD.close()
        query_string=''.join(["%ID=@",noLD_out]) #will query the file with the no LD buddy SNPs
        p=subprocess.Popen(["bcftools","query","-f","%ID:%CHROM:%POS\n","-i",query_string,vcf],shell=False,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        so,se=p.communicate()
        coord_list=so.splitlines()

        #To Do: how many variants don't have LD buddies and are also not the VCF that we check for coordinates
        # as it stands now we will skip these when querying because we don't have coordinates 
        if len(error_list)!=len(coord_list): #check that these match up
            print >> sys.stderr, "rsIDs for markers without LD buddies are not equal to coordinates found\n"

        with open(fn, 'a') as output:
            for coord in coord_list:
                rsid,chrom,pos=coord.split(":")
                output.write("\t".join([rsid,"".join(["chr",chrom,":",pos]),"\t".join(["NA"]*8),rsid]))
                output.write("\n")

        output.close()

    return 0

        
def get_proxy(token,pop,minrsq,rsid):
    """ make list of lists for the query for variants meeting minimum rsq"""
    try:
        r=requests.get(''.join(["https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=",rsid,"&pop=",pop,"&r2_d=r2&token=",token]))
    except requests.exceptions.RequestException as e:
        print >> sys.stderr, "Exception found %s" % e
        return 20
#    print(r.text    query_list=r.text.split("\n")
    query_list=r.text.split("\n")
    #print(query_list)
    for check in query_list: #handle exceptions when request is successful but rsID is not present
        if 'error' in check:
            print >> sys.stderr, "Error found in request output. rsID is probably not in 1000G reference"
            return 20
    

    #parse query results
    new_query_list=[]
    for entry in query_list:
        entry_list=entry.split("\t")
        if entry_list[0]!='': #check for empty last line
            if len(new_query_list)==0:
                new_query_list.append(entry_list) #keep header
            elif float(entry_list[6]) >= minrsq: #keep results meeting min rsq 
                new_query_list.append(entry_list)
    return new_query_list #list of all variants meeting minimum rsq

    #use Brooke's API token for  https://ldlink.nci.nih.gov/?tab=apiaccess
    #curl -k -X GET 'https://ldlink.nci.nih.gov/LDlinkRest/ldproxy?var=rs3&pop=MXL&r2_d=r2&token=5281fa13d3e3'
#does it need coordinattes or rsids

def check_data(bcf,out,vcf):
    """ Function to check what SNPs of interest and their proxy SNPs exist in HUNT data"""
    fn=".".join([out,"proxy.txt"]) #the out file where the results were written from read_rsID
    command=open_zip(fn)
    count=0
    marker_bed = NamedTemporaryFile(delete=True,suffix=".bed") #make temporary file for query
    marker_dict={}
    with open(marker_bed.name, 'w') as tmp:
        with command as f:
            for line in f:
                ls=line.rstrip()
                line_list=ls.split("\t")
                if count==0: #skip header
                    count=count+1
                    next
                else:
                    coords=line_list[1].split(":")
                    chrom=coords[0].replace("chr","")
                    tmp.write("\t".join([chrom,str(int(coords[1])-1),coords[1]])) #write bed file with coordinate
                    marker_dict[line_list[1]]=[line_list[10],line_list[6]] #save to dictionary to reference later 
        f.close()
        
    #bcftools query
    bcf_fn=".".join([out,"bcftools.txt"])
    subprocess.call([bcf,"query",vcf,"-R",marker_bed.name,"-f","%CHROM\t%POS\t%ID\n","-o",bcf_fn])
    tmp.close()

    print(marker_dict)
    
    #match back up with coords so we know what snp of interest it pertains to
    final_fn=".".join([out,"study.txt"])
    with open(final_fn,'w+') as final: #open file to write 
        bcommand=open_zip(bcf_fn)
        with bcommand as bcf_output: #open bcf query  output 
            for line in bcf_output:
                ls=line.rstrip()
                line_list=ls.split("\t")
                coordinate="".join(["chr",line_list[0],":",line_list[1]])
                if coordinate in marker_dict.keys():
                    print(marker_dict[coordinate])
                    snp_of_interest=marker_dict[coordinate][0] #match up with snp of interest rsID
                    r2=marker_dict[coordinate][1]
                    final.write("\t".join(["\t".join(line_list),snp_of_interest,r2]))
                    final.write("\n")
                
    bcf_output.close()
    final.close()
                
    # subprocess.call([bcftools,"query",vcf,"-R",tmp.name,"-f","%ID\t%CHROM\t%POS\t%REF\t%ALT[\t%DS]\n","-o",outName])
    return 0

###############################
######## MAIN ###############
############################
            
if __name__ == '__main__':
    args = argparser.parse_args()
    print >> sys.stderr, args
    
    if args.in_rs is not None and args.in_coord is not None:
        raise Exception('Supply rsID or coordinates but not both\n')
    elif args.in_rs is not None:
       read_rsid(args.in_rs,args.dbsnp_vcf,args.convert,args.token,args.population,args.minrsq,args.output) #rsID input
    elif args.in_coord is not None:
        read_coord(args.in_rs,args.dbsnp_vcf,args.convert) #coordinate input 
    else:
        raise Exception('Supply either rsID or coordinate file\n')

    #check HUNT data for snps of interest and their proxy snps, output a list of SNP IDs which can be used with a downstream script to pull out genotypes 
    check_data(args.bcftools,args.output,args.study_vcf)

 #TO DO: add open targets API 
