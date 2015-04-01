[![Travis-CI Build Status](https://travis-ci.org/thibautjombart/apex.png?branch=master)](https://travis-ci.org/thibautjombart/apex)


apex: Phylogenetic Methods for Multiple Gene Data
=================================================

The following functions enable the export from *apex* to other packages:
* **multidna2genind:** concatenates genes and export to genind

This is illustrated below:

```r
## find source files in apex
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)

## import data
x <- read.multiFASTA(files)
x
```

```
## === multidna ===
## [ 24 DNA sequences in 4 genes ]
## 
## @n.ind: 8 individuals
## @n.seq: 24 sequences in total
## @labels: 2340_50156.ab1  2340_50149.ab1  2340_50674.ab1  2370_45312.ab1  2340_50406.ab1  2370_45424.ab1 ...
## 
## @dna:
## $patr_poat43
## 8 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 764 
## 
## Labels: 2340_50156.ab1  2340_50149.ab1  2340_50674.ab1  2370_45312.ab1  2340_50406.ab1  2370_45424.ab1  ...
## 
## Base composition:
##     a     c     g     t 
## 0.320 0.158 0.166 0.356 
## 
## $patr_poat47
## 8 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 626 
## 
## Labels: 2340_50156.ab1  2340_50149.ab1  2340_50674.ab1  2370_45312.ab1  2340_50406.ab1  2370_45424.ab1  ...
## 
## Base composition:
##     a     c     g     t 
## 0.227 0.252 0.256 0.266 
## 
## $patr_poat48
## 8 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 560 
## 
## Labels: 2340_50156.ab1  2340_50149.ab1  2340_50674.ab1  2370_45312.ab1  2340_50406.ab1  2370_45424.ab1  ...
## 
## Base composition:
##     a     c     g     t 
## 0.305 0.185 0.182 0.327 
## 
## $patr_poat49
## 8 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 556 
## 
## Labels: 2340_50156.ab1  2340_50149.ab1  2340_50674.ab1  2370_45312.ab1  2340_50406.ab1  2370_45424.ab1  ...
## 
## Base composition:
##     a     c     g     t 
## 0.344 0.149 0.187 0.320
```

```r
## export to genind
obj <- multidna2genind(x)
obj
```

```
## 
##    #####################
##    ### Genind object ### 
##    #####################
## - genotypes of individuals - 
## 
## S4 class:  genind
## @call: DNAbin2genind(x = concatenate(x, genes = genes))
## 
## @tab:  8 x 22 matrix of genotypes
## 
## @ind.names: vector of  8 individual names
## @loc.names: vector of  11 locus names
## @loc.nall: number of alleles per locus
## @loc.fac: locus factor for the  22 columns of @tab
## @all.names: list of  11 components yielding allele names for each locus
## @ploidy:  1 1 1 1 1 1
## @type:  codom
## 
## Optional contents: 
## @strata: - empty -
## @hierarchy:  - empty -
## @pop:  - empty -
## @pop.names:  - empty -
## 
## @other: - empty -
```

