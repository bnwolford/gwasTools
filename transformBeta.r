#!/usr/bin/Rscript

#This script converts betas from linear mixed models to beta/se/OR comparable to logistic regression output for use in meta-analysis

### Lloyd Jones et al Genetics 2018
#http://www.genetics.org/content/genetics/early/2018/02/02/genetics.117.300360.full.pdf
#http://cnsgenomics.com/shiny/LMOR/
source("shiny_lmor_func.R")
library(data.table)
file<-"CAD.intermediate.all.imputed.results.toTransform.txt"
dt<-fread(file)
df<-data.frame(dt)
cases<-7710
controls<-58577
prev<-cases/(cases+controls)

#use standard error transformation (OR2 in paper)
res<-LmToOddsRatio(df,prev,1)
res$BETA_SE<-log(res$OR_SE) #convert OR to Beta
output<-"CAD.intermediate.all.imputed.results.toTransform_Vischer_transformed_SE.txt"
write.table(res, file=output,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

#use allele frequencies for transformation (OR1 in paper)
res<-LmToOddsRatio(df,prev,0)
res$BETA_AF<-log(res$OR) 
output<-"CAD.intermediate.all.imputed.results.toTransform_Vischer_transformed_AF.txt"
write.table(res, file=output,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)

#### Cook et al EJHG 2017
#system('python2 convert_beta_se.py', wait=FALSE)

