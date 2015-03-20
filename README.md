[![Travis-CI Build Status](https://travis-ci.org/thibautjombart/apex.png?branch=master)](https://travis-ci.org/thibautjombart/apex)

#apex
Extension of the R package ape for multiple genes.

Using the *apex* package
-------------
To install the development version from github: 
```r
library(devtools)
install_github("thibautjombart/apex")
```

The main (formal/S4) class is *multidna*.
See ?multidna for more information.
Here is a short example:
```r
library("apex")
## empty object
new("multidna")

## simple conversion with nicely ordered output
data(woodmouse)
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
x <- new("multidna", genes)
x
par(mfrow=c(3,1))
image(woodmouse)
image(x@dna[[1]])
image(x@dna[[2]])

## trickier conversion with missing sequences / wrong order
genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[c(5:1,14:15),501:965])
x <- new("multidna", genes)
x
par(mfrow=c(2,1))
image(x@dna[[1]])
image(x@dna[[2]])

```

