#!/usr/bin/Rscript

# Copyright (c) 2019 Brooke Wolford
# Revised from Dr. Lars Fritsche
# Lab of Dr. Cristen Willer and Dr. Mike Boehnke
# University of Michigan

options(stringsAsFactors=F)
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")
library("ggplot2")
library("tidyr")

option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited, can be gzipped, requires MAF and PVALUE columns"),   
  make_option("--prefix", type="character", default="",
    help="Prefix of output files"),   
  make_option("--top.size", type="numeric", default=0.125,
    help="top size = proportion of total length y axis [default=0.125]"),
  make_option("--break.top", type="numeric", default=15,
    help="set axis break at -log10(P) [default=15]"),
  make_option("--width", type="numeric", default=900,
    help="Width QQ plot in pixel [default=900]"),
  make_option("--height", type="numeric", default=900,
    help="Height QQ plot in pixel [default=900]"),
  make_option("--pointsize", type="numeric", default=16,
              help="Point size of plots [default=16]"),
  make_option("--pvalue", type="character", default="Pvalue",
              help="name of column with p.value [default='Pvalue']"),
  make_option("--gene",type="character",default="Gene",
              help="name of column with gene name [default='Gene']"),  
  make_option("--log10p", type="logical", default=F,
    help="Input p.value column with -log10(p.value) [default=F]"),    
  make_option("--maintitle", type="character", default="",
              help="Plot title"),
  make_option("--significance",type="numeric",
              help="Significance threshold. Bonferroni on tests used if not provided."),
  make_option("--pdf",type="logical",default=F,
    help="Plot as pdf [default=F]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list, description="This script creates qqplots for SAIGE-GENE output\n")

#################################################
################ FUNCTIONS ######################
#################################################


# QQ plot function
qqplotdata <- function(logpvector){
    o = sort(logpvector,decreasing=T)
    e = -log10(ppoints(length(o)))
    qqdata <- data.frame(o,e)
    qqdata$o <- round(qqdata$o,3)
    qqdata$e <- round(qqdata$e,3)
    keepU <- which(!duplicated(qqdata))
    qqdata <- qqdata[keepU,]
    
    N <- length(logpvector) ## number of p-values
    ## create the confidence intervals
    qqdata$c975 <- NA
    qqdata$c025 <- NA

            ## the jth order statistic from a
            ## uniform(0,1) sample
            ## has a beta(j,n-j+1) distribution
            ## (Casella & Berger, 2002,
            ## 2nd edition, pg 230, Duxbury)

    for(i in 1:length(keepU)){
        j <- keepU[i]
        qqdata$c975[i] <- -log10(qbeta(0.975,j,N-j+1))
        qqdata$c025[i] <- -log10(qbeta(0.025,j,N-j+1))
    }
    return(qqdata)
}

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
lambdaGC<-function(log10P){
    denom<-qchisq(0.5, df=1) #calculate denominator
    char<-sapply(log10P,log10toP) #convert from log10P to character(P) vector
    numer<-sapply(char,function(x) {as.numeric(x)}) #convert to numeric vector
    #print(summary(numer)) #print summary of p-values
    num<-qchisq(median(numer),df=1,lower.tail=F) #calculate numerator
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

#horizontal lines and corresponding colors
colLine <- c("red")

#open file, even if zipped
if (grepl('.gz',opt$input)) {
    gwas <- fread(paste(sep=" ","zcat",opt$input),header=T)
} else {
    gwas <- fread(opt$input, header=T)
}

## write out top 50 genes
order<-gwas[order(gwas[[opt$pvalue]]),]
signif<-head(data.frame(order[,c(opt$gene,opt$pvalue),with=F]),50)
signif[[opt$pvalue]]<-as.numeric(signif[[opt$pvalue]])
write.table(format(signif,digits=3,scientific=TRUE),file=paste(sep=".",opt$prefix,"top50.txt"),sep="\t",quote=FALSE,row.name=FALSE,col.name=TRUE)

#convert pvalues to -log10pvalue or use existing values in that scale
if(!opt$log10p) {
    gwas[[opt$pvalue]]<- as.numeric(gwas[[opt$pvalue]]) #handle any NAs
    gwas$log10P <- -log10(gwas[[opt$pvalue]])
    ycol <- "log10P"
} else { 
    ycol <- opt$pvalue
}

gwas<-gwas[complete.cases(gwas),] #remove NAs

ntest<-nrow(gwas) #number of genes tested
print(ntest)
if (is.null(opt$significance)){
    ##bonferroni line on number of tests 
    bon<-0.05/ntest
    yLine<-c(-log10(bon))
} else {
    yLine<-c(-log10(opt$significance))
}



## histogram of number of markers per category
cat_col<-grep("MACCate",names(gwas))
if (length(cat_col) > 1){
    d<-data.frame(apply(data.frame(gwas[,cat_col,with=F]), 2, as.numeric))
    d$sum<-rowSums(d)


    if (opt$pdf==TRUE) { #plot as pdf, default for height/width/point size are customized for png
        pdf(filename = paste0(opt$prefix,"_hist.pdf"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
    } else {
        png(filename = paste0(opt$prefix,"_hist.png"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
    }
    print(ggplot(gather(d), aes(value)) + 
        geom_histogram(bins = 30) +  facet_wrap(~key, scales = 'free_x')) + theme_bw()
    dev.off()

    }



#subset to gene name and p.value
gwas <- na.omit(data.frame(gwas[,c(opt$gene,ycol),with=F]))


## Generate QQ plot data by frequency bin
fbin <- character(0)
fN <- integer(0)
fx <- numeric(0)
fy <- numeric(0)


plotdata <- qqplotdata(gwas[[ycol]])
lambda<-lambdaGC(gwas[[ycol]]) #calculate lambda for this bin
print(format(lambda,digits=3))

fx <- c(fx,plotdata$e)
fy <- c(fy,plotdata$o)
conf <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                        'y'=c(plotdata$c975,rev(plotdata$c025)))
color<-"#56B4E9"

## QQ plot by binned frequencies
if (opt$pdf==TRUE) { #plot as pdf, default for height/width/point size are customized for png
    pdf(filename = paste0(opt$prefix,"_QQ.pdf"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
} else {
    png(filename = paste0(opt$prefix,"_QQ.png"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
}
	xlim <- c(0,max(fx,na.rm=T))
	ylim <- c(0,max(fy,na.rm=T))
	maxY <- max(fy,na.rm=T)
	par(mar=c(5.1,5.1,4.1,1.1))
	# plot version with two axes
	if(maxY > opt$break.top){
		# create pretty y-axis labels
        lab1 <- pretty(c(0,opt$break.top),n=ceiling(12 * (1-opt$top.size)))
        lab1 <- c(lab1[lab1 < opt$break.top],opt$break.top)
        lab2 <- pretty(c(opt$break.top,maxY),n=max(3,floor(12 * opt$top.size)))
        lab2 <- lab2[lab2 > max(lab1)]

        # resulting range of top scale in bottom scale units
        top.range = opt$break.top/(1 - opt$top.size) - opt$break.top
        top.data = max(lab2)-opt$break.top
        
        # function to rescale the top part
            rescale = function(y) { opt$break.top+(y-opt$break.top)/(top.data/top.range)}
            rescaled.y = rescale(fy[fy>opt$break.top])
            plot(0,0,
                 ylim=c(min(fy),opt$break.top*(1+opt$top.size)),xlim=xlim,axes=FALSE,
                 xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
                 ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
                 cex=1,cex.lab=1.5,cex.axis=1.5,bty="n",col="transparent",
                 main=opt$maintitle,pch=19)
		
		# Plot confidence intervals	

            polygon(conf$'x',ifelse(conf$'y'>opt$break.top,rescale(conf$'y'),conf$'y'),
                    col=grDevices::rgb(t(grDevices::col2rgb(color)),alpha=50,maxColorValue=255),border = NA)

		# add points
            points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)

            ## To do: add gene labels
            
		# identify line & add axis break
		lines(xlim,xlim,col="black",lty = 2)
		axis(1,cex.axis=1.5,cex.lab=1.5)
		par(las=1)
		axis(side=2,at=lab1,cex.axis=1.5,cex.lab=1.5)
		par(las=0)
		box()
		par(las=0)
		points(fx[fy>opt$break.top],rescaled.y,cex=1,col=fcol[fy>opt$break.top],pch=19)
		par(las=1)
		axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=1.5,cex.lab=1.5)
		axis.break(axis=2,breakpos=opt$break.top,style="zigzag",brw=0.02)
		axis.break(axis=4,breakpos=opt$break.top,style="zigzag",brw=0.02)
		lines(range(fx),c(opt$break.top,opt$break.top),col = "grey",lty = 6)
		abline(h=ifelse(yLine<opt$break.top,
			yLine,
			rescale(yLine)),
			col=colLine,lwd=1.5,lty=2)

	# plot version with single y axes
	} else {
		par(mar=c(5.1,5.1,4.1,1.1),las=1)
		axislim <- ceiling(range(xlim,ylim,yLine))
		plot(0,0,
			ylim=axislim,xlim=xlim,axes=T,
			xlab=expression(plain(Expected)~~group("(",-log[10]*italic(P),")")),
			ylab=expression(plain(Observed)~~group("(",-log[10]*italic(P),")")),
			cex=1,cex.lab=1.5,cex.axis=1.5,col="transparent",
			main=opt$maintitle,pch=19)
		# Plot confidence intervals

                polygon(conf$'x',conf$'y',
					col=grDevices::rgb(t(grDevices::col2rgb(color)),alpha=50,maxColorValue=255),
					border = NA)
                ##add points 
		points(fx,fy,col=color,pch=19)
		# identity line & genome-wide significance line
		lines(axislim,axislim,col = "grey",lwd=1.5,lty=2)
		abline(h=yLine,col=colLine,lwd=1.5,lty=2)

                text(0,0,paste(sep=" ",expression(lambda),format(lambda,digits=3)))
                
	}
dev.off()


