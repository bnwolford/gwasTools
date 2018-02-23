options(stringsAsFactors=F)
library("plotrix")
library("data.table")
library("RColorBrewer")
library("optparse")

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
    qqdata$c95 <- NA
    qqdata$c05 <- NA

    ## the jth order statistic from a
    ## uniform(0,1) sample
    ## has a beta(j,n-j+1) distribution
    ## (Casella & Berger, 2002,
    ## 2nd edition, pg 230, Duxbury)

    for(i in 1:length(keepU)){
        j <- keepU[i]
        qqdata$c95[i] <- -log10(qbeta(0.95,j,N-j+1))
        qqdata$c05[i] <- -log10(qbeta(0.05,j,N-j+1))
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
    print(summary(numer)) #print summary of p-values
    num<-qchisq(median(numer),df=1,lower.tail=F) #calculate numerator
    lam<-num/denom #calculate lambda
    return(lam)
}
    
option_list <- list(
  make_option("--input", type="character", default="",
    help="Input file, tab delimited"),   
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
  make_option("--maintitle", type="character", default="",
    help="Plot title"),
  make_option("--pdf",type="logical",default=F,
    help="Plot as pdf [default=F]")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

# horizontal lines and corresponding colors
yLine <- c(-log10(5E-8))
colLine <- c("red")

#read file
gwas <- fread(opt$input)

#convert pvalues to -log10pvalue or use existing values in that scale
if(!opt$log10p) {
    gwas[[opt$pvalue]]<- as.numeric(gwas[[opt$pvalue]]) #handle any NAs
    gwas$log10P <- -log10(gwas[[opt$pvalue]])
    ycol <- "log10P"
} else { 
    ycol <- opt$pvalue
}

gwas<-gwas[complete.cases(gwas),] #remove NAs

#establish maf column for qqplot
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
    }
} else if (opt$minMAC > 0) { #minMAC provided so filter by MAC
    if (opt$mac %in% colnames(gwas)) {
        gwas<-gwas[gwas[[opt$mac]] > opt$minMAC] #filter by MAC
    } else if (opt$ac %in% colnames(gwas)) {
        gwas$mac<-gwas[[opt$ac]]
        gwas$mac[which(gwas$mac > gwas[[opt$sample.size]])] <- 2*gwas[[opt$sample.size]][which(gwas$mac > gwas[[opt$sample.size]])] - gwas$mac[which(gwas$mac > gwas[[opt$sample.size]])] #convert to MAC
        gwas<-gwas[mac > opt$minMAC] #filter by MAC
    } else {
        stop("Please provide --ac or --mac arguments that match a column in input file\n")
    }
} else {
    warning("Results are not filtered by MAF or MAC\n")
}

# TO DO: plot  bins by MAC rather than MAF

#subset to maf and p.value
gwas <- na.omit(data.frame(gwas[,c(mafcol,ycol),with=F]))

# Determine frequency bins and create variable for binned QQ plot
freqbins <- c(c(0.5,0.05,0.005,0.001,0)[which(c(0.5,0.05,0.005,0.001,0) > floor(minMAF*1000000)/1000000)],floor(minMAF*1000000)/1000000)
gwas$freqbin <- cut(gwas[[mafcol]], freqbins, include.lowest=T)
freqtable <- table(gwas$freqbin)
freqtable <- freqtable[order(-as.numeric(gsub("[\\[\\(](.+),.+","\\1",names(freqtable))))]
freqtable <- freqtable[freqtable > 0]

#initialize variables to return lambda values per MAF bin
lambda_file_name<-paste0(opt$prefix,"_lambda.txt")
lambda_df<-NULL

## Generate QQ plot data by frequency bin
fbin <- character(0)
fN <- integer(0)
fx <- numeric(0)
fy <- numeric(0)
fcol <- character(0)
legendcol <- character(0)
conf <- list()
allcols <- brewer.pal(4,"Set1")
for(f in 1:length(freqtable)){
	fbin <- c(fbin,names(freqtable)[f])
	fsnps <- which(gwas$freqbin ==names(freqtable)[f])
	plotdata <- qqplotdata(gwas[[ycol]][fsnps])
        lambda<-lambdaGC(gwas[[ycol]][fsnps]) #calculate lambda for this bin
        lambda_df<-rbind(lambda_df,data.frame(lambda=lambda,frequency_bin=fbin[f])) #make lambda data frame
	fN <- c(fN,freqtable[f])
	fx <- c(fx,plotdata$e)
	fy <- c(fy,plotdata$o)
	fcol <- c(fcol,rep(allcols[f],length(plotdata$o)))
	conf[[f]] <- data.frame('x'=c(plotdata$e,rev(plotdata$e)),
                                'y'=c(plotdata$c95,rev(plotdata$c05)))
	legendcol <- c(legendcol,allcols[f])
}
legendtext <- paste0("MAF=",fbin,"; N SNPs=",format(fN,big.mark=",",scientific=FALSE))

#write data frame of lambda values
write.table(x=lambda_df,file=lambda_file_name,col.names=T,row.names=F,quote=F,sep="\t")

## QQ plot by binned frequencies
if (opt$pdf==TRUE) { #plot as pdf, default for height/width/point size are customized for png
    pdf(filename = paste0(opt$prefix,"_QQ.png"), width = opt$width, height = opt$height, pointsize = opt$pointsize)
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
		for(p in 1:length(conf)){
			polygon(conf[[p]]$'x',ifelse(conf[[p]]$'y'>opt$break.top,rescale(conf[[p]]$'y'),conf[[p]]$'y'),
				col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
				border = NA)
		}

		# add points
		points(fx[fy<opt$break.top],fy[fy<opt$break.top],cex=1,col=fcol[fy<opt$break.top],pch=19)

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
		legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
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
		for(p in 1:length(conf)){
				polygon(conf[[p]]$'x',conf[[p]]$'y',
					col=grDevices::rgb(t(grDevices::col2rgb(allcols[p])),alpha=50,maxColorValue=255),
					border = NA)
		}
		points(fx,fy,col=fcol,pch=19)
		# identity line & genome-wide significance line
		lines(axislim,axislim,col = "grey",lwd=1.5,lty=2)
		abline(h=yLine,col=colLine,lwd=1.5,lty=2)
		legend("topleft",legend=legendtext,col=legendcol,pch=15,bty="n")
	}
dev.off()
