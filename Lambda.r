#!/usr/bin/Rscript

# Copyright (c) 2018 Brooke Wolford
# Revised from Dr. Lars Fritsche
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")


option_list <- list(
  make_option("--input", type="character", default="",
              help="Input file, tab delimited, can be gzipped, requires MAF and PVALUE columns"),   
  make_option("--prefix", type="character", default="",
              help="Prefix of output files"),   
  make_option("--maf", type="character", default='MAF',
              help="name of column with MAF [default='MAF']"), 
  make_option("--af",type="character", default='AF',
              help="name of column with AF [default='AF']"),
  make_option("--mac",type="character",default="MAC",
              help="name of column with MAC [default='MAC']"),
  make_option("--ac",type="character",default="AC",
              help="name of column with AC [default='AC']"),
  make_option("--sample.size",type="character",default="N",
              help="name of column with sample size, required to convert allele count to to minor allele count [default='N']"),
  make_option("--minMAF",type="numeric", default=0,
              help="minimum MAF of variants for plotting [default=0]"),
  make_option("--minMAC",type="numeric",default=0,
              help="minimum MAC of variants for plotting [default=0]"),
  make_option("--pvalue", type="character", default="PVALUE",
              help="name of column with p.value [default='PVALUE']"),
  make_option("--log10p", type="logical", default=F,
              help="Input p.value column with -log10(p.value) [default=F]"),
  make_option("--exclude",type="character",default="",
              help="File with SNP IDs to exclude. Assumes input file has SNP ID")
              
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script creates qqplots across MAF bins for summary statistics and calculates lambda genomic control at percentiles of the chi squared distribution.\n")

#################################################
################ FUNCTIONS ######################
#################################################

# convert -log10(P) values to as.character(P)
log10toP <- function(log10P){
    log10P <- abs(as.numeric(log10P))
    if(is.na(log10P)) return(NA)
    if(log10P==Inf) return(as.character(0))
    if(log10P > 300){
        part1 <- log10P%/%100*100
        part2 <- log10P-part1
        P <- format(signif(10^-part2,6), scientific = T)
        P <- paste(as.numeric(gsub("e-.+","",P)),"E-",as.numeric(gsub(".+-","",P),sep="")+part1,sep="")
    } else {
        P <- signif(10^-log10P,6)
    }
    return(as.character(P))
}


#calculate lambda for genomic correction
lambdaGC<-function(log10P,qval){
    denom<-qchisq(qval, df=1) #calculate denominator
    char<-sapply(log10P,log10toP) #convert from log10P to character(P) vector
    numer<-sapply(char,function(x) {as.numeric(x)}) #convert to numeric vector
    num<-qchisq(quantile(numer,qval),df=1,lower.tail=F) #calculate numerator
    lam<-num/denom #calculate lambda
    return(lam)
}

################################################
############## MAIN #############################
#################################################

#parse arguments
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#TO DO: check for required inputs

#open file, even if zipped
if (grepl('.gz',opt$input)) {
    gwas <- fread(paste(sep=" ","zcat",opt$input),header=T)
} else {
    gwas <- fread(opt$input, header=T)
}

#convert pvalues to -log10pvalue or use existing values in that scale
if(!opt$log10p) {
    gwas[[opt$pvalue]]<- as.numeric(gwas[[opt$pvalue]]) #handle any NAs
    gwas$log10P <- -log10(gwas[[opt$pvalue]])
    ycol <- "log10P"
} else { 
    ycol <- opt$pvalue
}

gwas<-gwas[complete.cases(gwas),] #remove NAs

#establish maf column
if (opt$maf %in% colnames(gwas)) { 
    mafcol<-opt$maf
    gwas[[opt$maf]]<-as.numeric(gwas[[opt$maf]])
    minMAF<-min(gwas[[opt$maf]]) #minMAF of all input data
} else if (opt$af %in% colnames(gwas)) {
    gwas$maf<-as.numeric(gwas[[opt$af]])
    gwas$maf[which(gwas$maf > 0.5)] <- 1 - gwas$maf[which(gwas$maf > 0.5)] #convert AF to MAF
    mafcol<-"maf"
    minMAF<-min(gwas[[mafcol]]) #minMAF of all input data
} else {
    stop("Please provide --af or --maf arguments that match a column in input file\n")
}

#filter by minMAF or minMAC if provided
if (opt$minMAF > 0) { #minMAF provided so filter by MAF
    if (opt$minMAC > 0 ) {
        stop("Please only provide either --minMAF or --minMAC but not both\n")
    } else {
        gwas <- gwas[gwas[[mafcol]] > opt$minMAF] #filter by MAF
        minMAF<-min(gwas[[mafcol]]) #new minMAF
        print(summary(gwas[[mafcol]])) #print MAF 
    }
} else if (opt$minMAC > 0) { #minMAC provided so filter by MAC
    if (opt$mac %in% colnames(gwas)) {
        gwas<-gwas[gwas[[opt$mac]] > opt$minMAC] #filter by MAC
        print(summary(gwas[[opt$mac]])) #print MAC
    } else if (opt$ac %in% colnames(gwas)) {
        gwas$mac<-gwas[[opt$ac]]
        gwas$mac[which(gwas$mac > gwas[[opt$sample.size]])] <- 2*gwas[[opt$sample.size]][which(gwas$mac > gwas[[opt$sample.size]])] - gwas$mac[which(gwas$mac > gwas[[opt$sample.size]])] #convert to MAC
        gwas<-gwas[mac > opt$minMAC] #filter by MAC
        print(summary(gwas$mac)) #print MAC
    } else {
        stop("Please provide --ac or --mac arguments that match a column in input file\n")
    }
} else {
    warning("Results are not filtered by MAF or MAC\n")
}

## TO DO:  bins by MAC rather than MAF

## To Do : use exclude file to remove SNPs with true signal

##subset to maf and p.value
gwas <- na.omit(data.frame(gwas[,c(mafcol,ycol),with=F]))

## Determine frequency bins and create variable for binned QQ plot
freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
gwas$freqbin <- cut(gwas[[mafcol]], freqbins, include.lowest=T)
freqtable <- table(gwas$freqbin)
freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
freqtable <- freqtable[freqtable > 0]

##initialize variables to return lambda values per MAF bin and per percentile
lambda_file_name<-paste0(opt$prefix,"_percentiles_lambda.txt")
lambda_df<-NULL


fbin <- character(0)

for(f in 1:length(freqtable)){
    for (q in seq(.1,.9,0.1)){
        fbin <- c(fbin,names(freqtable)[f])
        fsnps <- which(gwas$freqbin ==names(freqtable)[f])
        lambda<-lambdaGC(gwas[[ycol]][fsnps],q) #calculate lambda for this bin
        lambda_df<-rbind(lambda_df,data.frame(lambda=lambda,frequency_bin=fbin[f],quantile=q)) #make lambda data frame
}}

##write data frame of lambda values
write.table(x=lambda_df,file=lambda_file_name,col.names=T,row.names=F,quote=F,sep="\t")

