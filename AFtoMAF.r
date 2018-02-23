options(stringsAsFactors=F)
library(data.table)
library(optparse)

#this script converts the AF to MAF and outputs a file <prefix>_minor.txt with a new column

option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited"),   
  make_option("--prefix", type="character", default="",
    help="Prefix of output files"),   
  make_option("--af",type="character",default="AF",
    help="name of column with AF [default='AF']"),
  make_option("--colName",type="character",default="MAF",
    help="name of new column with MAF [default='MAF']") 
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#read file
file <- fread(opt$input)

if (opt$af %in% colnames(file)) { #check maf column exists
    file$maf<-as.numeric(file[[opt$af]]) #make new column
    file$maf[which(file$maf > 0.5)] <- 1 - file$maf[which(file$maf > 0.5)] #convert AF to MAF
} else {
    stop("Please provide --af argument that match a column in input file\n")
}

colnames(file)[colnames(file)=="maf"] <- opt$colName

#write file
filename=paste0(opt$prefix,"_minor.txt")
write.table(x=file,file=filename,col.names=T,row.names=F,quote=F,sep="\t")

