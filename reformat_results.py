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
parser.add_argument("-p","--pheno",help="Phenotye file to identify cases and controls from VCF",type=str,required=True)
parser.add_argument("-v","--vcf",help="VCF for a specific chromosome",type=str,required=True)
parser.add_argument("-c","--chr",help="Chromosome of the VCF, 23 for X",type=int,required=True)
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
            info_dict[coord]={}
            info_dict[coord]["hwe"]=lineList[6]
            info_dict[coord]["rsq"]=lineList[3]
    return(info_dict)

# read phenotype file and get sample IDs that are cases or controls, output as 2 lists
def read_pheno(pheno_file):
    case=[]
    control=[]
    #read cases and controls from phenotype file
    openCommand=open_file(pheno_file)
    with openCommand as f:
        for line in f:
            ls=line.rstrip()
            lineList=ls.split("\t")
            if lineList[12] == "1": #cases
                case.append(lineList[0])
            elif lineList[12] == "0": #controls
                control.append(lineList[0])
            else:
                print >> sys.stderr, lineList[0] + "has no phenotype\n"
    return case,control

# get relevant per variant information from the VCF
def parse_VCF(vcf,case,control):
    #get sample list
    proc=subprocess.Popen(["/usr/local/bin/bcftools", "query","-l",vcf],stdout=subprocess.PIPE)
    query=filter(lambda x:len(x)>0,(line.strip() for line in proc.stdout)) #get list of samples
    
    count_dict={}
    openCommand=open_file(vcf)
    with openCommand as f:
        for line in f:
            ls=line.rstrip()
            if not ls.startswith('#'): #skip headers
                lineList=ls.split("\t")
                coord=":".join([lineList[0],lineList[1]])
                genos=lineList[9:]
                vcf_dict=dict(zip(query,genos)) #per variant dictionary of samples and genotypes
                out=geno_counts(vcf_dict,case,control) #make counts for each variant 
                count_dict[coord]=out #save list of counts to the count dictionary on a per variant basis
    return count_dict

#turn 0|0:0.000 into counts for cases and controls
def geno_counts(dict,case,control):
    ncase=len(case)
    ncontrol=len(control)
    n0_control=0 
    n1_control=0
    n2_control=0
    n0_case=0
    n1_case=0
    n2_case=0
    for sample in dict.keys(): #evaluate every sample's genotype in terms of effect allele
        geno=dict[sample].split(":")[0]
        if sample in control:
            if geno=="1|1":
                n2_control+=1
            elif geno=="0|0":
                n0_control+=1
            elif geno=="0|1" or geno=="1|0":
                n1_control+=1
            #else:
                #print >> sys.stderr, "Geno %s not recognized.\n" %geno
        elif sample in case:
            if geno=="1|1":
                n2_case+=1
            elif geno=="0|0":
                n0_case+=1
            elif geno=="0|1" or geno=="1|0":
                n1_case+=1
            #else:
                #print >> sys.stderr, "Geno %s not recognized.\n" %geno
        #else:
            #print >> sys.stderr, "Sample %s not a case or control.\n" %sample
    n0=n0_control+n0_case #homozygous non effect
    n1=n1_control+n1_case #heterozygous
    n2=n2_control+n2_case #homozygous effect
    eaf=(n2*2 + n1) / ((ncase+ncontrol)*2)
    eaf_case=(n2_case*2 + n1_case)/(ncase*2)
    eaf_control=(n2_control*2 + n1_control)/(ncontrol*2)
    return([n0,n1,n2,eaf,eaf_case,eaf_control,n0_control,n1_control,n2_control,n0_case,n1_case,n2_case])

def reformat_file(file,info_dict,vcf_dict,chrom):
    openCommand=open_file(file)
    with openCommand as f:
        count=0
        for line in f: #iterate over every line
            if count==0:
                header=line
                count+=1
                new_header="\t".join(["SNP","STRAND","BUILD","CHR","POS","EFFECT_ALLELE","NON_EFFECT_ALLELE","N","N0","N1","N2","EAF","N0_controls","N1_controls","N2_controls","EAF_controls","N0_cases","N1_cases","N2_cases","EAF_cases","HWE","CALL_RATE","BETA","SE","PVAL","IMPUTED","INFO"])
                print(new_header)
                count+=1
            else:
                ls = line.rstrip()
                lineList=ls.split("\t")
                if str(lineList[1]) == str(chrom): #restrict to per chrom so we can parallelize and use a single chrom VCF for info
                    coord=":".join([lineList[1],lineList[2]])
                    #pull info from other sources
                    strand="+"
                    build="37"
                    call_rate="1.000"
                    imputed="1"
                    effect_allele=lineList[4]
                    non_effect_allele=lineList[3]
                    try:
                        hwe=info_dict[coord]["hwe"]
                        rsq=info_dict[coord]["rsq"]
                        n0,n1,n2,eaf,eaf_case,eaf_control,n0_control,n1_control,n2_control,n0_case,n1_case,n2_case=vcf_dict[coord] #pull info out of vcf
                    except KeyError:
                        hwe="."
                        rsq="."
                        n0,n1,n2,eaf,eaf_case,eaf_control,n0_control,n1_control,n2_control,n0_case,n1_case,n2_case=["."]*12
                        pass
                    new_line="\t".join([lineList[0],strand,build,lineList[1],lineList[2],effect_allele,non_effect_allele,str(lineList[7]), n0, n1, n2, lineList[6],n0_control, n1_control, n2_control,eaf_control, n0_case, n1_case, n2_case, eaf_case,hwe,call_rate,lineList[8],lineList[9],lineList[11],imputed,rsq])
                    print(new_line)
                    count+=1
            

#########################
########## MAIN #########
#########################

def main():
                                   
    info_dict=info(args.info) # read in hwe and rsq from a separate file
    case,control=read_pheno(args.pheno)
    vcf_dict=parse_VCF(args.vcf,case,control) #parse n and allele frequency for all, case, control
    reformat_file(args.file,info_dict,vcf_dict, args.chr)
    
if __name__ == "__main__":
    main()
        
