#!/usr/bin/Rscript

# Copyright (c) 2018 Brooke Wolford
# Revised from Dr. Sarah Graham
# Implements published methods, see below
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

option_list <- list(
    make_option("--input", type="character", default=NULL,
                help="Input file, tab delimited, can be gzipped"),
    make_option("--k",type="character",
                help="Name of column with sample prevalance. Otherwise prevalence from --numCase and --numControl will be used for all genetic variants."),
    make_option("--beta",type="character",default="BETA",
                help="Name of column with beta [default=BETA]"),
    make_option("--se",type="character",default="SE",
                help="Name of column with se [default=SE]"),
    make_option("--freq",type="character",default="FREQ",
                help=" Name of column with freq [default=FREQ]"),
    make_option("--n",type="character",default="N",
                help="Name of column with sample size N [default=N]"),
    make_option("--stdErrTrans",type="logical",default=TRUE,
                help="Logical for performing Lloyd-Jones et al standard error transformation (requires SE, N, BETA) [default=TRUE"),
    make_option("--afTrans",type="logical",default=TRUE,
                help="Logical for performing Lloyd-Jones et al allele frequency based transformation (requires BETA, FREQ, and k or numCases and numControls) [default=TRUE]"),
    make_option("--cook",type="logical",default=TRUE,
                help="Logical for performing Cook et al transformation (requires BETA, SE, numCase, numControl) [default=TRUE]"),
    make_option("--pathToFunc",type="character",default="shiny_lmor_func.R",
                help="Full path to shiny_lmor_func.R from Lloyd-Jones et al"),
    make_option("--numCase", type="numeric", default=NULL,
                help="Number of cases"),
    make_option("--numControl",type="numeric", default=NULL,
                help="Number of controls"),
    make_option("--prefix",type="character",default=NULL,
                help="Output file prefix")
      )
parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script transforms summary statistics from linear mixed models to be compatible with those from logistic/linear regression for the purposes of meta-analysis.\n")

##################################################################
########################## FUNCTIONS #############################
##################################################################

#### Cook et al EJHG 2017
#doi:10.1038/ejhg.2016.150
#Guidance for the utility of linear models in meta-analysis of genetic association studies of binary phenotypes
cook_function<-function(cases,controls,beta,se,out_prefix){
    cat("Performing transformation based on Cook et al doi:10.1038/ejhg.2016.150\n")
    outCook<-paste0(out_prefix,"_Cook.txt")
    frac_cases<-cases/(cases+controls)
    frac_controls<-controls/(cases+controls)
    df$BETA_T<-df[[beta]]/(frac_cases * frac_controls) #transformed beta
    df$SE_T<-df[[se]]/(frac_cases * frac_controls) #transformed SE
    write.table(df,file=outCook, sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}


### Lloyd Jones et al Genetics 2018
#doi:10.1534/genetics.117.300360
#Transformation of Summary Statistics from Linear Mixed Model Association on All-or-None Traits to Odds Ratio
#http://cnsgenomics.com/shiny/LMOR/ hosts shiny_lmor_func.R
#derive OR based on beta and standard error
lj_se<-function(se,n,beta,prev,out_prefix){
    cat("Performing transformation based on Lloyd-Jones et al doi:10.1534/genetics.117.300360\n")
    source(opt$pathToFunc)
    outLJ_SE<-paste0(out_prefix,"_LloydJones_SE.txt")

    #rename data columns to something that shiny_lmor_func.R expects
    names(df)[names(df) == se] <- "SE"
    names(df)[names(df) == beta] <- "BETA"
    names(df)[names(df) == n] <- "N"

    res<-LmToOddsRatio(df,prev,1) #execute function
    res$BETA_SE<-log(res$OR_SE) #convert OR to Beta
    write.table(res, file=outLJ_SE,sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

#derive OR based on beta and allele frequencies
lj_af<-function(beta,freq,prev,out_prefix){
    cat("Performing transformation based on Lloyd-Jones et al doi:10.1534/genetics.117.300360\n")
    source(opt$pathToFunc)
    outLJ_AF<-paste0(out_prefix,"_LloydJones_AF.txt")

    #rename data columns to something that shiny_lmor_func.R expects
    names(df)[names(df) == freq] <- "FREQ"
    names(df)[names(df) == beta] <- "BETA"
    
    res<-LmToOddsRatio(df,prev,0) #execute function
    res$BETA_AF<-log(res$OR) #convert OR to Beta
    write.table(res, file=outLJ_AF, sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)
}

#function to check if argument exists
exists<-function(arg){
    if (sum(opt == as.factor(arg)) > 0) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#######################################################
##################### MAIN #############################
#######################################################


#parse arguments
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#check for input file and output file arguments
if (!exists(opt$input) || !exists(opt$prefix)){
    stop("Please provide input file --input and output file --output\n")
}

#read in file, even if gzipped
file<-opt$input
if (grepl('.gz',file)){
    cat("Reading gzipped file\n.")
    dt<-fread(paste(sep=" ","zcat",file),header=T)
} else {
    cat("Reading file.\n")
    dt<-fread(file) #read file
}
df<-data.frame(dt) #convert to data frame which LmToOddsRatio expects


#### perform Cook et al transformation
if (opt$cook) {
    #check for required arguments
    if (!exists(opt$numCase) || !exists(opt$numControl) || !exists(opt$beta) || !exists(opt$se)) {
        stop("Please provide number of cases and controls, name of beta column and name of se column for Cook et al transformation\n")
    } else {
        cook_function(opt$numCase,opt$numControl,opt$beta,opt$se,opt$prefix)
    }
}

#### perform Llyod Jones et al standard error transformation
if (opt$stdErrTrans) {

    #check for arguments required for prevalence
    if (exists(opt$k)){
        prev<-df[[opt$k]] #use column with prevalence per variant
    } else if (exists(opt$numCase) && exists(opt$numControl)) {
        prev<-opt$numCase/(opt$numCase+opt$numControl) #use prevalence for all genetic variants
    } else {
        stop("Please provide number of cases and controls or name of column with prevalence for Lloyd-Jones et al transformation\n")
    }
    
    #check for other required arguments
    if (!exists(opt$se) || !exists(opt$n) || !exists(opt$beta)) {
        stop("Please provide name of columns with SE, N, and BETA for Lloyd-Jones et al standard error transformation\n")
    } else {
        lj_se(opt$se, opt$n, opt$beta, prev, opt$prefix)
    }
}

####perform Lloyd-Jones et al allele frequency transformation
if (opt$afTrans) {

    #check for arguments required for prevalence
    if (exists(opt$k)){
        prev<-df[[opt$k]] #use column with prevalence per variant
    } else if (exists(opt$numCase) && exists(opt$numControl)) {
        prev<-opt$numCase/(opt$numCase+opt$numControl) #use prevalence for all genetic variants
    } else {
        stop("Please provide number of cases and controls or name of column with prevalence for Lloyd-Jones et al transformation\n")
    }

    #check for other required arguments
    if (!exists(opt$beta) || !exists(opt$freq) ) {
        stop("Please provide name of columns with SE and BETA for Lloyd-Jones et al allele frequency based transformation\n")
    } else {
        lj_af(opt$beta,optfreq,prev,opt$prefix)
    }
}








