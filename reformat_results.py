#!/usr/bin/env python

#input file column order
"""1SNPID
2CHR
3POS
4Allele0
5Allele1
6AC
7AF
8N
9BETA
10SE
11Tstat
12p.value
13p.value.NA
14Is.SPA.converge
15varT
16varTsta
"""
############################
##### IMPORT MODULES #######
###########################

import os
import sys
import subprocess
from subprocess import check_call
import argparse
import gzip

###########################
##### PARSE ARGUMENTS ####
##########################
parser = argparse.ArgumentParser(description="File to format sex-straifeid analysis results in the format requested\n")
parser.add_argument("-f","--file",help="File of SAIGE results to reformat",type=str, required=True)
parser.add_argument("-i","--info",help="File with information like hwe and rsq",type=str,required=True)
parser.add_argument("-ka","--plink_case",help=".hwe file for cases",type=str,required=True)
parser.add_argument("-ko","--plink_control",help=".hwe file for controls",type=str,required=True)
parser.add_argument("-c","--chr",help="Chromosome of the VCF, 23 for X",type=int,required=True)
parser.add_argument("-d","--header",help="Is there a header on -f file",action='store_true')
args=parser.parse_args()

###############################
######## FUNCTIONS  #########
###############################

def open_file(file):
    if ".gz" in file:
        command=gzip.open(file,"rt")
        print >> sys.stderr, "Opening gzipped/bgzipped file %s" %file
    elif file == "-":
        command=sys.stdin
        print >> sys.stderr, "Opening file %s from stdin" %file
    elif ".gz" not in file:
        command=open(file,"rt")
        print >> sys.stderr, "Opening file %s" %file
    return command

# read info file with columns chr, posE, SNPID, rsq, MAF,SNPID, hwe p-value
def info(info_file):
    info_dict={} #initalize dictionary for hwe, info score
    openCommand=open_file(info_file)
    with openCommand as f:
        for line in f:
            ls=line.rstrip()
            lineList=ls.split("\t")
            coord=":".join([lineList[0],lineList[1]])
            snp=lineList[2] #use snp id instead of coord to account for multiallelic variants
            info_dict[snp]={}
            info_dict[snp]["hwe"]=lineList[6]
            info_dict[snp]["rsq"]=lineList[3]
    return(info_dict)

#read in .hwe files for cases and controls from plink
def geno_counts(case,control):
    geno_dict={}
    open_case=open_file(case)
    with open_case as a:
        next(a)
        for line in a:
            ls=line.rstrip()
            lineList=ls.split()
            geno_dict[lineList[1]]=lineList[3:6] #use SNP ID instead of coord for unique alleles
    open_control=open_file(control)
    with open_control as o:
        next(o)
        for line in o:
            ls=line.rstrip()
            lineList=ls.split()
            geno_dict[lineList[1]].append(lineList[5])
    return(geno_dict)

def parse_counts(list,effect_allele):
    minor,major,case,control=list
    if major==effect_allele: #alt is effect allele
        n0_case,n1_case,n2_case=case.split("/")
        n0_control,n1_control,n2_control=control.split("/") #plink hwe are hom minor, het, hom major
    elif minor==effect_allele:
        n2_case,n1_case,n0_case=case.split("/")
        n2_control,n1_control,n0_control=control.split("/")
    else:
        print >> sys.stderr, "Effect allele is neither minor or major allele from PLINK files\n"
    n0_case,n1_case,n2_case,n0_control,n1_control,n2_control=map(int,[n0_case,n1_case,n2_case,n0_control,n1_control,n2_control]) #convert strings to integers 
    n0=n0_case+n0_control
    n1=n1_case+n1_control
    n2=n2_case+n2_control
    ncase=n0_case+n1_case+n2_case
    ncontrol=n0_control+n1_control+n2_control
    eaf_case=(n2_case*2 + n1_case)/(ncase*2)
    eaf_control=(n2_control*2 + n1_control)/(ncontrol*2)
    #n0,n1,n2,eaf_case,eaf_control,n0_control,n1_control,n2_control,n0_case,n1_case,n2_case
    parsed_list=map(str,[n0,n1,n2,eaf_case,eaf_control,n0_control,n1_control,n2_control,n0_case,n1_case,n2_case])
    return(parsed_list) #return as strings

def reformat_file(file,info_dict,geno_dict,chrom,header):
    new_header="\t".join(["SNP","STRAND","BUILD","CHR","POS","EFFECT_ALLELE","NON_EFFECT_ALLELE","N","N0", "N1","N2","EAF","N0_controls","N1_controls","N2_controls","EAF_controls","N0_cases","N1_cases","N2_cases","EAF_cases","HWE","CALL_RATE","BETA","SE","PVAL","IMPUTED","INFO"])
    print(new_header)
    
    openCommand=open_file(file)
    with openCommand as f:
        if header==True: #skip first line if there is a header
            next(f)
        for line in f: #iterate over every line
            ls = line.rstrip()
            lineList=ls.split() #any white space delimiter
            if str(lineList[1]) == str(chrom): #restrict to per chrom so we can parallelize and use a single chrom VCF for info
                coord=":".join([lineList[1],lineList[2]])
                snp=lineList[0]
                #pull info from other sources
                strand="+"
                build="37"
                call_rate="1.000"
                imputed="1"
                effect_allele=lineList[4]
                non_effect_allele=lineList[3]
                eaf=lineList[6]
                n=lineList[7]
                try:
                    hwe=info_dict[snp]["hwe"]
                    rsq=info_dict[snp]["rsq"]
                    n0,n1,n2,eaf_case,eaf_control,n0_control,n1_control,n2_control,n0_case,n1_case,n2_case=parse_counts(geno_dict[snp],effect_allele)
                except KeyError:
                    hwe="."
                    rsq="."
                    n0,n1,n2,eaf_case,eaf_control,n0_control,n1_control,n2_control,n0_case,n1_case,n2_case=["."]*11
                    pass
                new_line="\t".join([lineList[0],strand,build,lineList[1],lineList[2],effect_allele,non_effect_allele,str(n), n0, n1, n2, str(eaf),n0_control, n1_control, n2_control,eaf_control, n0_case, n1_case, n2_case, eaf_case,hwe,call_rate,lineList[8],lineList[9],lineList[11],imputed,rsq])
                print(new_line)
            

#########################
########## MAIN #########
#########################

def main():
                                   
    info_dict=info(args.info) # read in hwe and rsq from a separate file
    geno_dict=geno_counts(args.plink_case,args.plink_control)
    reformat_file(args.file,info_dict,geno_dict, args.chr,args.header)
    
if __name__ == "__main__":
    main()
        
