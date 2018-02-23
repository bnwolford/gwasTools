options(stringsAsFactors=F)
library(data.table)
library(optparse)

option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited"),   
  make_option("--prefix", type="character", default="",
    help="Prefix of output files"),   
  make_option("--ac",type="character",default="AC",
    help="name of column with AC [default='AC']"),
  make_option("--sample.size",type="character",default="N",
    help="name of column with sample size, required to convert allele count to to minor allele count [default='N']")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

#raed file
file <- fread(opt$input)

#establish mac column
if (opt$ac %in% colnames(file)) {
    file$mac<-file[[opt$ac]]
    file$mac[which(file$mac > file[[opt$sample.size]])] <- 2*file[[opt$sample.size]][which(file$mac > file[[opt$sample.size]])] - file$mac[which(file$mac > file[[opt$sample.size]])] #convert to MAC
} else {
    stop("Please provide --ac argument that match a column in input file\n")
}

#write file
filename=paste0(opt$prefix,"_minor.txt")
write.table(x=file,file=filename,col.names=T,row.names=F,quote=F,sep="\t")

