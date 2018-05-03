#!/usr/bin/Rscript

#Copyright (c) 2018 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
library(data.table)
library(optparse)

option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited, can be gzipped"),   
  make_option("--prefix", type="character", default="",
    help="Output file name prefix"),   
  make_option("--pvalue",type="character",default="MAC",
    help="name of column with pvalue [default='p.value']"),
  make_option("--log10p", type="logical", default=F,
    help="Input p.value column with -log10(p.value) [default=F]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script prints the min, 1st quartile, median, mean, 3rd quartile, and max of the negative log 10 pvalues from a summary statistics file\n")

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#check for required arguments
if (opt$input=="" || opt$prefix=="") {
    stop("Please provide --input and --output arguments\n")
}

#open file, even if zipped
if (grepl('.gz',opt$input)) {
    file <- fread(paste(sep=" ","zcat",opt$input),header=T)
} else {
    file <- fread(opt$input, header=T)
}

col<-opt$pvalue
library(data.table)

options(digits=8) #set digits for summary

if(!opt$log10p){
    file[[col]]<-as.numeric(file[[col]]) #handle NA
    file$neglog10<--log10(file[[col]]) #convert to neg log P
    summ<-summary(file[file$neglog10!=Inf,]$neglog10) #summary
} else{
    summ<-summary(file[[col]]) #summary since the value are already neglog10
}

#write output
filename<-paste(sep="_",opt$prefix,"pval_summary.txt")
cat("P-value summary", summ, file=filename,sep="\n",append=FALSE)


