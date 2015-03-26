[![Travis-CI Build Status](https://travis-ci.org/thibautjombart/apex.png?branch=master)](https://travis-ci.org/thibautjombart/apex)





#apex
Extension of the R package ape for multiple genes.

Installing *apex*
-------------
To install the development version from github: 

```r
library(devtools)
install_github("thibautjombart/apex")
```

Functionalities
----------------

#### Classes of object

See the classes:
* **multidna:** formal (S4) class, storing data using a list of DNAbin objects.
* 
Example code:

```r
library("apex")
```

```
## Error in library("apex"): there is no package called 'apex'
```

```r
## empty object
new("multidna")
```

```
## Error in getClass(Class, where = topenv(parent.frame())): "multidna" is not a defined class
```

```r
## simple conversion with nicely ordered output
data(woodmouse)
```

```
## Warning in data(woodmouse): data set 'woodmouse' not found
```

```r
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
```

```
## Error in eval(expr, envir, enclos): object 'woodmouse' not found
```

```r
x <- new("multidna", genes)
```

```
## Error in getClass(Class, where = topenv(parent.frame())): "multidna" is not a defined class
```

```r
x
```

```
## Error in eval(expr, envir, enclos): object 'x' not found
```

```r
par(mfrow=c(3,1))
image(woodmouse)
```

```
## Error in image(woodmouse): object 'woodmouse' not found
```

```r
image(x@dna[[1]])
```

```
## Error in image(x@dna[[1]]): object 'x' not found
```

```r
image(x@dna[[2]])
```

```
## Error in image(x@dna[[2]]): object 'x' not found
```

```r
## trickier conversion with missing sequences / wrong order
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[c(5:1,14:15),501:965])
```

```
## Error in eval(expr, envir, enclos): object 'woodmouse' not found
```

```r
x <- new("multidna", genes)
```

```
## Error in getClass(Class, where = topenv(parent.frame())): "multidna" is not a defined class
```

```r
x
```

```
## Error in eval(expr, envir, enclos): object 'x' not found
```

```r
plot(x)
```

```
## Error in plot(x): object 'x' not found
```

#### Reading data from multiple files
See the functions:
* **read.multidna:** reads multiple DNA alignments with various formats
* **read.multiFASTA:** same for FASTA files

Example code:

```r
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
files
```

```
## character(0)
```

```r
## read files
x <- read.multiFATSA(files)
```

```
## Error in eval(expr, envir, enclos): could not find function "read.multiFATSA"
```

```r
x
```

```
## Error in eval(expr, envir, enclos): object 'x' not found
```

```r
plot(x)
```

```
## Error in plot(x): object 'x' not found
```



#### Data handling
See the functions:
* **concatenate:** concatenate seeral genes into a single DNAbin matrix 
* **x[i,j]:** subset x by individuals (i) and/or genes (j)

Example code:

```r
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
files
```

```
## character(0)
```

```r
## read files
x <- read.multiFASTA(files)
```

```
## Error in eval(expr, envir, enclos): could not find function "read.multiFASTA"
```

```r
x
```

```
## Error in eval(expr, envir, enclos): object 'x' not found
```

```r
plot(x)
```

```
## Error in plot(x): object 'x' not found
```

```r
## subset
plot(x[1:3,2:4])
```

```
## Error in plot(x[1:3, 2:4]): object 'x' not found
```

```r
## concatenate
y <- concatenate(x)
```

```
## Error in eval(expr, envir, enclos): could not find function "concatenate"
```

```r
y
```

```
## Error in eval(expr, envir, enclos): object 'y' not found
```

```r
image(y)
```

```
## Error in image(y): object 'y' not found
```


#### Exporting data
Check functions:
* **multidna2genind:** concatenate genes and export to genind

Example:

```r
## read data in
files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
files
```

```
## character(0)
```

```r
## read files
x <- read.multiFASTA(files)
```

```
## Error in eval(expr, envir, enclos): could not find function "read.multiFASTA"
```

```r
## export to genind
obj <- multidna2genind(x)
```

```
## Error in eval(expr, envir, enclos): could not find function "multidna2genind"
```

```r
obj
```

```
## Error in eval(expr, envir, enclos): object 'obj' not found
```
