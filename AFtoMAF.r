#!/usr/bin/Rscript

#Copyright (c) 2018 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
library(data.table)
library(optparse)

#this script converts the AF to MAF and outputs a file <prefix>_minor.txt with a new column

option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited, can be gzipped"),   
  make_option("--output", type="character", default="",
    help="Name for output file"),   
  make_option("--af",type="character",default="AF",
    help="name of column with AF [default='AF']"),
  make_option("--colName",type="character",default="MAF",
    help="name of new column with MAF [default='MAF']") 
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script converts allele frequency to minor allele frequency in a new column titled MAF or --colName and writes the output to a new file ")

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

#calcualte maf from af    
if (opt$af %in% colnames(file)) { #check maf column exists
    file$maf<-as.numeric(file[[opt$af]]) #make new column
    file$maf[which(file$maf > 0.5)] <- 1 - file$maf[which(file$maf > 0.5)] #convert AF to MAF
} else {
    stop("Please provide --af argument that match a column in input file\n")
}

#rename column to MAF or provided --colName
colnames(file)[colnames(file)=="maf"] <- opt$colName

#write file
filename<-opt$output
write.table(x=file,file=filename,col.names=T,row.names=F,quote=F,sep="\t")
