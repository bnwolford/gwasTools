# gwasTools

A collection of R scripts that might be useful for exploring and plotting GWAS results.
Derived from scripts from Dr. Lars Fritsche, Dr. Wei Zhou, Dr. Sarah Graham.

The following R packages need to be installed for running these Rscripts:

*optparse, data.table, RColorBrewer, plotrix, pwr, ggplot2*


## QQ Plot

Frequency binned and thinned QQ plot. Also calculates median genomic control inflation factor, lambda. Please check out the required/available Rscript parameters by using the following command:

    Rscript QQplot.r --help

<img src="/Example_QQ.png" height="480">

## Manhattan Plot

Fast Manhattan plot script. Please check out the required/available Rscript parameters by using the following command:

    Rscript ManhattanPlot.r --help

<img src="/Example_Manhattan.png" height="480">


## Minimal/example Input format for QQ and Manhattan plots

|CHROM	|POS	|MAF	|PVALUE	|
|---	|---	|---	|---	|
|1  	|1  	|0.05	|0.99	|
|2  	|2   	|0.15	|0.1	|
|3  	|3  	|0.5	|0.25	|

## ACtoMAC.r 

Convert alelle count (AC) to minor allele count (MAC) in additional column titled MAC or custom name (--colName). Requires columns with AC and sample size (N) and an input and output file name.  
   
## AFtoMAF.r 

Convert allele frequency (AF) to minor allele frequency (MAF) in additional column titled MAF or custom name (--colName). Requires columns with AF and an input and output file name.
    
## comparisonPlot.R 

Create scatterplots to compare values from two columns in a data frame. Print Pearson's correlation of the value sto standard out. Please check out the required/available Rscript parameters by using the following command:

    Rscript comparisonPlot.R --help
    
<img src="/Example_comparison.png" height="480">
    
## transformBeta.r 

Transform beta values from linear mixed model so they are comparable to beta values from logistic/linear regression. This implements work from Cook et al (doi:10.1038/ejhg.2016.150) and Lloyd-Jones et al (doi:10.1534/genetics.117.300360). Please check out the required/available Rscript parameters by using the following command:

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
    
    Rscript AFtoMAF.r
    --input ExampleGWAS.txt \
    --output Example_minor.txt \
    --af A1Freq \
    --colName MAF
 
    Rscript ACtoMAC.r
    --input ExampleGWAS.txt \
    --output Example_minor.txt \
    --ac AC \
    --colName MAC

    Rscript comparisonPlot.R
    --input ExampleGWAS.txt \
    --prefix Example_comparison \
    --VAL1 p.value
    --VAL2 Logreg.pval
    --labdata1 "SAIGE" \
    --labdata2 "HAIL" \
    --title "Pvalue comparison in all variants" \
    --negLog10 TRUE

    Rscript transformBeta.R
    --input ExampleGWAS.txt \ 
    --prefix ExampleGWAS \
    --numCase 1000 \
    --numControl 4000 \
    --pathToFunc shiny_lmor_func.R
