
##
## Internal functions: not documented, not exported
##
## Thibaut Jombart, April 2015
##

## compute the number of missing sequences for each gene in a list of DNAbin matrices 'x'
.nMissingSequences <- function(x){
    ##if(!inherits(x, "multidna")) stop("x is not a multidna object")
    out <- sapply(x, function(e) sum(apply(e==as.DNAbin("-"),1,all)))
    return(sum(out,na.rm=TRUE))
}
