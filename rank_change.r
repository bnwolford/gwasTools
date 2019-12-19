#!/usr/bin/Rscript

#Copyright (c) 2018 Brooke Wolford
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
library(data.table)
library(optparse)
library(ggplot2)
library(ggrepel)

option_list <- list(
  make_option("--file1", type="character", default="",
              help="Input file1, tab delimited, can be gzipped"),
  make_option("--file2",type="character",default="",
              help="Input file2, tab delimited, can be gzipped"),
  make_option("--prefix", type="character", default="",
              help="Output file name prefix"),   
  make_option("--pvalue1",type="character",default="p.value",
              help="name of column with pvalue in file 1[default='p.value']"),
  make_option("--pvalue2",type="character",default="p.value",
              help="name of column with pvalue in file 2[default='p.value']"),
  make_option("--merge",type="character",default="Gene",
              help="name of column to merge file1 and file2 on[default='Gene']")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script prints the min, 1st quartile, median, mean, 3rd quartile, and max of the negative log 10 pvalues from a summary statistics file\n")

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#check for required arguments
if (opt$file1=="" || opt$file2=="" ||  opt$prefix=="") {
    stop("Please provide --file1 and --file2 and --prefix arguments\n")
}

##function to open file even if zipped
open<-function(file){
    if (grepl('.gz',file)) {
        file <- fread(paste(sep=" ","zcat",file),header=T)
    } else {
        file <- fread(file, header=T)
    }
    return(file)
}

col1<-opt$pvalue1
col2<-opt$pvalue2
file1<-open(opt$file1)
file2<-open(opt$file2)

#print(head(file1))
#print(head(file2))

options(digits=8) #set digits for summary

file1<-file1[order(file1[[col1]])]
file2<-file2[order(file2[[col2]])]
file1$rank<-seq.int(nrow(file1))
file2$rank<-seq.int(nrow(file2))

df<-merge(file1,file2,by=opt$merge)
df$diff<-ifelse(df$rank.x-df$rank.y>1000,TRUE,FALSE)
#sub<-head(df,100)



#write output
filename<-paste(sep="_",opt$prefix,"rankChange.pdf")
pdf(filename,height=5,width=5,useDingbats=FALSE)
ggplot(df,aes(x=rank.x,y=rank.y)) + geom_point(alpha=0.5) + theme_bw() + geom_text_repel(data=df[df$diff==TRUE,],aes(label=get(opt$merge)),size=2) + geom_abline(linetype="dashed",color="red",slope=1,intercept=0) + labs(x=opt$file1,y=opt$file2)
dev.off()


