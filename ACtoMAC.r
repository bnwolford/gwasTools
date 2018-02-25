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
  make_option("--output", type="character", default="",
    help="Prefix of output files"),   
  make_option("--ac",type="character",default="AC",
    help="name of column with AC [default='AC']"),
  make_option("--sample.size",type="character",default="N",
    help="name of column with sample size, required to convert allele count to to minor allele count [default='N']"),
  make_option("--colName",type="character",default="MAC",
    help="name of new column with MAC [default='MAC']")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script converts allele count to minor allele count in a new column titled MAC or --colName, and writes the output to a new file.\n")

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#check for required arguments
if (opt$input=="" || opt$output=="") {
        stop("Please provide --input and --output arguments\n")
}

#open file, even if zipped
if (grepl('.gz',opt$input)) {
    file <- fread(paste(sep=" ","zcat",opt$input),header=T)
} else {
    file <- fread(opt$input, header=T)
}

#establish mac column
if (opt$ac %in% colnames(file)) {
    file$mac<-file[[opt$ac]]
    file$mac[which(file$mac > file[[opt$sample.size]])] <- 2*file[[opt$sample.size]][which(file$mac > file[[opt$sample.size]])] - file$mac[which(file$mac > file[[opt$sample.size]])] #convert to MAC
} else {
    stop("Please provide --ac argument that matches a column in input file\n")
}

#rename column to MAC or provided --colName
colnames(file)[colnames(file)=="mac"]<-opt$colName

#write file
filename<-opt$output
write.table(x=file,file=filename,col.names=T,row.names=F,quote=F,sep="\t")

