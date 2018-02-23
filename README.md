# gwasTools

A collection of R scripts that might be useful for exploring and plotting GWAS results.

The following R pakages need to be installed for running these Rscripts:

*optparse, data.table, RColorBrewer, plotrix, pwr, ggplot2*


## QQ Plot

Frequency binned and thinned QQ plot.

Please check out the required/available Rscript parameters by using the following command

    Rscript QQplot.r --help

<img src="/Example_QQ.png" height="480">

## Manhattan Plot

Fast Manhattan plot script. Please check out the required/available Rscript parameters by using the following command

    Rscript ManhattanPlot.r --help

<img src="/Example_Manhattan.png" height="480">


## Minimal/example Input format for QQ and Manhattan plots

|CHROM	|POS	|MAF	|PVALUE	|
|---	|---	|---	|---	|
|1  	|1  	|0.05	|0.99	|
|2  	|2   	|0.15	|0.1	|
|3  	|3  	|0.5	|0.25	|

## ACtoMAC.r 

    Rscript ACtoMAC.r --help
    
## AFtoMAF.r 

    Rscript AFtoMAF.r --help
    
## comparisonPlot.R 

    Rscript comparisonPlot.R --help
    
## transformBeta.r 

    Rscript transformBeta.r --help 
    


### Examples (run from command line):

    Rscript QQplot.r \
    --input ExampleGWAS.txt \
    --prefix Example \
    --maf MAF \
    --pvalue PVALUE \
    --maintitle 'An Example QQ plot'


    Rscript ManhattanPlot.r 
    --input ExampleGWAS.txt \
    --prefix Example \
    --chr CHROM \
    --pos POS \
    --pvalue PVALUE \
    --maintitle 'An Example Manhattan plot'
 
