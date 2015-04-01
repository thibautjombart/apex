[![Travis-CI Build Status](https://travis-ci.org/thibautjombart/apex.png?branch=master)](https://travis-ci.org/thibautjombart/apex)




*apex*: Phylogenetic Methods for Multiple Gene Data
=================================================
*apex* implements some new classes to handle DNA sequences from different genes and individuals.
It implements new classes extending object classes from *ape* and *phangorn* to store multiple gene data, and some useful wrappers mimicking existing functionalities for multiple genes.
This document provides an overview of the package's content.


Installing *apex*
-------------
To install the development version from github:

```r
library(devtools)
install_github("thibautjombart/apex")
```

The stable version can be installed from CRAN using:

```r
install.packages("apex")
```

Then, to load the package, use:

```r
library("apex")
```

```
## Loading required package: ape
## Loading required package: phangorn
```

```
## Warning: package 'phangorn' was built under R version 3.2.0
```


New object classes
------------------
Two new classes of object extend existing data structure for multiple genes:
* **multidna:** based on *ape*'s `DNAbin` class
* **xxxx:** based on ...

###  multidna
This formal (S4) class can be seen as a multi-gene extension of *ape*'s `DNAbin` class.
Data is stored as a list of DNAbin objects, with additional slots for extra information.
The class definition can be obtained by:
```r
getClassDef("multidna")
```
* **@dna**: list of `DNAbin` matrices, with matching rows
* **@labels**: labels of the individuals (rows of the matrices in `dna`)
* **@n.ind**: the number of individuals
* **@n.seq**: the total number of sequences in the dataset
* **@ind.info**: an optional dataset storing individual metadata 
* **@gene.info**: an optional dataset storing gene metadata 

Any of these slots can be accessed using `@` (see example below).

New `multidna` objects can be created via two ways:

1. using the constructor `new("multidna", ...)`
2. reading data from files (see section on 'importing data' below)

We illustrate the use of the constructor below (see `?new.multidna`) for more information.
We use *ape*'s dataset *woodmouse*, which we artificially split in two 'genes', keeping the first 500 nucleotides for the first gene, and using the rest as second gene. Note that the individuals need not match across different genes: matching is handled by the constructor.

```r
## empty object
new("multidna")
```

```
## === multidna ===
## [ 0 DNA sequence in 0 gene ]
## 
## @n.ind: 0 individual
## @n.seq: 0 sequence in total
## @labels:
```

```r
## using a list of genes as input
data(woodmouse)
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
x
```

```
## === multidna ===
## [ 30 DNA sequences in 2 genes ]
## 
## @n.ind: 15 individuals
## @n.seq: 30 sequences in total
## @labels: No305 No304 No306 No0906S No0908S No0909S...
## 
## @dna:
## $gene1
## 15 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 500 
## 
## Labels: No305 No304 No306 No0906S No0908S No0909S ...
## 
## Base composition:
##     a     c     g     t 
## 0.326 0.230 0.147 0.297 
## 
## $gene2
## 15 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 465 
## 
## Labels: No305 No304 No306 No0906S No0908S No0909S ...
## 
## Base composition:
##     a     c     g     t 
## 0.286 0.295 0.103 0.316
```

```r
## access the various slots
x@labels
```

```
##  [1] "No305"   "No304"   "No306"   "No0906S" "No0908S" "No0909S" "No0910S"
##  [8] "No0912S" "No0913S" "No1103S" "No1007S" "No1114S" "No1202S" "No1206S"
## [15] "No1208S"
```

```r
x@n.ind
```

```
## [1] 15
```

```r
class(x@dna) # this is a list
```

```
## [1] "list"
```

```r
names(x@dna) # names of the genes
```

```
## [1] "gene1" "gene2"
```

```r
x@dna[[1]] # first gene
```

```
## 15 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 500 
## 
## Labels: No305 No304 No306 No0906S No0908S No0909S ...
## 
## Base composition:
##     a     c     g     t 
## 0.326 0.230 0.147 0.297
```

```r
x@dna[[2]] # second gene
```

```
## 15 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 465 
## 
## Labels: No305 No304 No306 No0906S No0908S No0909S ...
## 
## Base composition:
##     a     c     g     t 
## 0.286 0.295 0.103 0.316
```

```r
## compare the input dataset and the new multidna
par(mfrow=c(3,1), mar=c(6,6,2,1))
image(woodmouse)
image(x@dna[[1]])
image(x@dna[[2]])
```

![plot of chunk class](vignettes/figs/class-1.png) 

```r
## same but with missing sequences and wrong order
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[c(5:1,14:15),501:965])
x <- new("multidna", genes)
x
```

```
## === multidna ===
## [ 22 DNA sequences in 2 genes ]
## 
## @n.ind: 15 individuals
## @n.seq: 22 sequences in total
## @labels: No305 No304 No306 No0906S No0908S No0909S...
## 
## @dna:
## $gene1
## 15 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 500 
## 
## Labels: No305 No304 No306 No0906S No0908S No0909S ...
## 
## Base composition:
##     a     c     g     t 
## 0.326 0.230 0.147 0.297 
## 
## $gene2
## 15 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 465 
## 
## Labels: No305 No304 No306 No0906S No0908S No0909S ...
## 
## Base composition:
##     a     c     g     t 
## 0.286 0.294 0.103 0.316
```

```r
par(mar=c(6,6,2,1))
plot(x)
```

![plot of chunk class](vignettes/figs/class-2.png) 



Importing data
--------------
Two simple functions permit to import data from multiple alignements into `multidna` objects:
* **read.multidna:** reads multiple DNA alignments with various formats
* **read.multiFASTA:** same for FASTA files

Both functions rely on the single-gene counterparts in *ape* and accept the same arguments.
Each file should contain data from a given gene, where sequences should be named after individual labels only.
Here is an example using a dataset from *apex*:

```r
## get address of the file within apex
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
files # this will change on your computer
```

```
## [1] "/usr/local/lib/R/site-library/apex/patr_poat43.fasta"
## [2] "/usr/local/lib/R/site-library/apex/patr_poat47.fasta"
## [3] "/usr/local/lib/R/site-library/apex/patr_poat48.fasta"
## [4] "/usr/local/lib/R/site-library/apex/patr_poat49.fasta"
```

```r
## read these files
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
names(x@dna) # names of the genes
```

```
## [1] "patr_poat43" "patr_poat47" "patr_poat48" "patr_poat49"
```

```r
par(mar=c(6,11,2,1))
plot(x)
```

![plot of chunk readfiles](vignettes/figs/readfiles-1.png) 



Handling data
--------------
Several functions facilitate data handling:
* **concatenate:** concatenate several genes into a single DNAbin matrix
* **x[i,j]:** subset x by individuals (i) and/or genes (j)

Example code:

```r
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
files
```

```
## [1] "/usr/local/lib/R/site-library/apex/patr_poat43.fasta"
## [2] "/usr/local/lib/R/site-library/apex/patr_poat47.fasta"
## [3] "/usr/local/lib/R/site-library/apex/patr_poat48.fasta"
## [4] "/usr/local/lib/R/site-library/apex/patr_poat49.fasta"
```

```r
## read files
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
par(mar=c(6,11,2,1))
plot(x)
```

![plot of chunk handling](vignettes/figs/handling-1.png) 

```r
## subset
plot(x[1:3,2:4])
```

![plot of chunk handling](vignettes/figs/handling-2.png) 

```r
## concatenate
y <- concatenate(x)
y
```

```
## 8 DNA sequences in binary format stored in a matrix.
## 
## All sequences of same length: 2506 
## 
## Labels: 2340_50156.ab1  2340_50149.ab1  2340_50674.ab1  2370_45312.ab1  2340_50406.ab1  2370_45424.ab1  ...
## 
## Base composition:
##     a     c     g     t 
## 0.298 0.187 0.197 0.319
```

```r
par(mar=c(5,8,2,1))
image(y)
```

![plot of chunk concat](vignettes/figs/concat-1.png) 


Exporting data
---------------
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

