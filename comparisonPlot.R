#!/usr/bin/Rscript

# Copyright (c) 2018 Brooke Wolford
# Revised from script by Dr. Wei Zhou
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
Sys.setlocale("LC_CTYPE", "C.UTF-8")

library(data.table)
library(optparse)
library(plotrix)
library(ggplot2)

option_list <- list(
    make_option(c("-p", "--prefix"), type="character", default="",
                help="prefix for output files"),
    make_option("--input", type="character", default="",
                help="input file, can be gzipped"),
    make_option("--VAL1", type="character", default="",
                help="header name for value 1 for comparison"),
    make_option("--VAL2", type="character", default="",
                help="header name for value 2 for comparison"),
    make_option("--labdata1", type="character", default="",
                help="data set name in xlab"),
    make_option("--labdata2", type="character", default="",
                help="data set name in ylab"),
    make_option("--title",type="character",default="",
                help="title"),
    make_option("--negLog10",type="logical", default=TRUE,
                help="Convert VAL1 and VAL2 to -log10 scale, will add -log10 to labels [default=TRUE]"),
    make_option(c("-c","--cor"),type="logical",default=TRUE,
                help="Calculate Pearson's correlation and print to stdout [default=TRUE]"),
    make_option(c("-f","--filter"),type="character",default="",
                help="Column name to filter on"),
    make_option("--max",type="numeric",default="",
                help="Maximum number to filter on (keep values <= this"),
    make_option("--min",type="numeric",default="",
                 help="Minimum number to filter on (keep values > to this")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script creates scatter plots comparing values from two columns of an input file (e.g. MAF vs Beta) and can print the Pearson's correlation to standard out.\n")
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
pheno <- opt$pheno
prefix <- opt$prefix
valcol1 <- opt$VAL1
valcol2 <- opt$VAL2
xlabdata <- opt$labdata1
ylabdata <- opt$labdata2
inputfile <- opt$input
title <- opt$title
log <- opt$negLog10
cor <- opt$cor

#check for required arguments
if (inputfile=="" || prefix =="" || valcol1=="" || valcol2=="") {
    stop("Please provide --input and --output arguments\n")
}

#open file, even if zipped
if (grepl('.gz',inputfile)) {
    data <- fread(paste(sep=" ","zcat",inputfile),header=T)
} else {
    data <- fread(inputfile, header=T)
}

#do we filter?
if (opt$filter != "") {
    print(nrow(data))
    if (opt$filter %in% colnames(data)) {
        if (opt$max != ""){
            data <- data[data[[opt$filter]] <= opt$max] #filter by max
            print(nrow(data))
        } else if (opt$min != ""){
            data <- data[data[[opt$filter]] > opt$min] #filter by min
            print(nrow(data))
        } else {
            stop("If you intend to filter on a column please provide --max or --min but not both\n")
        }
    }
    else {
        stop("Please provide --filter that is a column name of --input\n")
    }
}
    
#reset value names
setnames(data, c(valcol1, valcol2), c("VAL1", "VAL2"))

if (log == TRUE) {
    data$val1 <- -log10(data$VAL1)
    data$val2 <- -log10(data$VAL2)
    xlab<-paste0("-log10 ",xlabdata)
    ylab<-paste0("-log10 ",ylabdata)
} else {
    data$val1<-data$VAL1
    data$val2<-data$VAL2
    xlab<-xlabdata
    ylab<-ylabdata
}

# Pearson's correlation
if (cor==TRUE){
    correlation=cor(data[complete.cases(data),]$VAL1,data[complete.cases(data),]$VAL2)
    cat(correlation)
}

#plot
png(paste0(prefix, ".png"),width=800,height=800)
ggplot(data,aes(x=val1,y=val2)) + geom_point(alpha = 0.3) + theme_bw() + 
geom_abline(slope = 1, color="red",linetype="dashed", size=1.5,alpha=0.5) + 
theme(text = element_text(size=24),axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5)) +
labs(title=title,x=xlab,y=ylab)
dev.off()
