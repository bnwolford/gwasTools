#!/usr/bin/Rscript

#Copyright (c) 2019 Brooke Wolford
# Labs of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
library(data.table)
library(optparse)

option_list <- list(
    make_option("--old_bed", type="character", default=""),
    make_option("--new_bed",type="character",default=""),
    make_option("--unmapped",type="character",default=""),
    make_option("--prefix",type="character",default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script creates a key for old reference genome coords to new reference genome coordinates after liftover\n")

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

old<-fread(opt$old_bed)
new<-fread(opt$new_bed)
unmapped<-data.table(read.table(opt$unmapped,comment.char="#"))

names(old)<-c("old_chr","old_poss","old_pose")
names(new)<-c("new_chr","new_poss","new_pose")
names(unmapped)<-c("un_chr","un_poss","un_pose")

old$old_coord<-paste(sep=":",old$old_chr,old$old_poss)
new$new_coord<-paste(sep=":",new$new_chr,new$new_poss)
unmapped$coord<-paste(sep=":",unmapped$un_chr,unmapped$un_poss)

#initialize variables as we are adding columns to the old file
old$new_chr<-as.character(NA)
old$new_poss<-as.numeric(NA)
old$new_pose<-as.numeric(NA)
old$new_coord<-as.character(NA)

for (i in 1:nrow(old)){
    if (old[i]$old_coord %in% unmapped$coord){
        next
    } else{
        old[i]$new_chr<-new[i]$new_chr
        old[i]$new_poss<-new[i]$new_poss
        old[i]$new_pose<-new[i]$new_pose
        old[i]$new_coord<-new[i]$new_coord
    }
    
}
 
#write output
filename<-paste(sep=".",opt$prefix,"txt")
write.table(old, file=filename,sep="\t",append=FALSE,quote=FALSE)


