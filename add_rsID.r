#!/usr/bin/Rscript

#Copyright (c) 2018 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
library(data.table)
library(optparse)
library(tidyverse)

option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited, can be gzipped"),   
  make_option("--output", type="character", default="",
    help="Name for output file"),   
  make_option("--col",type="character",default="SNP",
    help="name of column with SNP ID formatted X:XXXXX_X/X [default='SNP']"),
  make_option("--dbsnp",type="character",default="",
    help="Bed file from dbsnp with columns chr, posS, posE, rsID")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script adds an rsID column to a .txt file with results from BOLT-LMM.")

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#check for required arguments
if (opt$input=="" || opt$output=="" || opt$dbsnp=="") {
    stop("Please provide --input and --output and --dbsnp arguments\n")
}

snp_col<-opt$col

#open file, even if zipped
if (grepl('.gz',opt$input)) {
    file <- fread(paste(sep=" ","zcat",opt$input),header=T)
} else {
    file <- fread(opt$input, header=T)
}
file<-as_tibble(file)

#open dbsnp file, even if zipped
if (grepl('.gz',opt$input)) {
    dbsnp <- fread(paste(sep=" ","zcat",opt$dbsnp),header=T)
} else {
    dbsnp <- fread(opt$dbsnp, header=T)
}
names(dbsnp)<-c("chr","posS","posE","rsID")
dbsnp<-as_tibble(dbsnp)

#split up SNP name in file
df_cols<-names(file)
#n_cols<-length(df_cols)
#snp_col<-which(df_cols==opt$col)
file<-separate(file, (!!snp_col), c("snp","alleles"),"_") %>% separate(snp,c("chr","posE"),":")
file<-mutate(file,chr=type.convert(chr)) %>% mutate(posE=type.convert(posE))

#inner join of file and dbsnp
join<-inner_join(file,dbsnp,by=c("chr"="chr","posE"="posE"))

#reformat join so we just added rsID column to the original data frame
join<-mutate(join,(!!snp_col):=paste(sep=":",chr,posE)) %>% mutate((!!snp_col):=paste(sep="_",SNP,alleles))
final<-select(join,one_of(c(df_cols,"rsID")))

#write file
filename<-opt$output
write.table(x=final,file=filename,col.names=T,row.names=F,quote=F,sep="\t")
