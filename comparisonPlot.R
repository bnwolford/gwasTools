#revised from Wei Zhou

options(stringsAsFactors=F)
Sys.setlocale("LC_CTYPE", "C.UTF-8")

library(data.table)
library(optparse)
library(plotrix)
library(ggplot2)
# Load extra R functions
source("/net/dumbo/home/zhowei/projects/UKBIOBANK/summary/script/helperFunctions.r")

option_list <- list(
    make_option(c("-p", "--prefix"), type="character", default="",
                help="prefix for output files"),
    make_option("--inputfile", type="character", default="",
                help="inputfile, can be gzipped"),
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
                help="Convert VAL1 and VAL2 to -log10 scale [default=TRUE]"),
    make_option(c("c","--cor"),type="logical",default=TRUE,
                help="Calculate Pearson's correlation and print to stdout [default=TRUE]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
pheno <- opt$pheno
prefix <- opt$prefix
valcol1 <- opt$VAL1
valcol2 <- opt$VAL2
xlabdata <- opt$labdata1
ylabdata <- opt$labdata2
inputfile <- opt$inputfile
title <- opt$title
log <- opt$negLog10
cor <- opt$cor

#open file, even if zipped
if (grepl('.gz',inputfile)) {
    data <- fread(paste(sep=" ","zcat",inputfile),header=T)
} else {
    data <- fread(inputfile, header=T)
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

#correlation
if (cor==TRUE){
    correlation=cor(data[complete.cases(data),]$VAL1,data[complete.cases(data),]$VAL2)
    print(correlation)
}

#plot
png(paste0(prefix, ".png"),width=800,height=800)
ggplot(data,aes(x=val1,y=val2)) + geom_point(alpha = 0.3) + theme_bw() + 
geom_abline(slope = 1, color="red",linetype="dashed", size=1.5,alpha=0.5) + 
theme(text = element_text(size=24),axis.text.x = element_text(size=24), axis.text.y = element_text(size=24), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"),plot.title = element_text(hjust = 0.5)) +
labs(title=title,x=xlab,y=ylab)
dev.off()
