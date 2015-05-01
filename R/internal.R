
##
## Internal functions: not documented, not exported
##
## Thibaut Jombart, April 2015
##

## compute the number of missing sequences for each gene in a list of DNAbin matrices or phyDat objects 'x'
.nMissingSequences <- function(x){
    ## only keep non-empty matrices
    x <- lapply(x, as.character)
    x <- x[sapply(x, nrow)>0]
    out <- sapply(x, function(e) sum(apply(e=="-",1,all)))
    return(sum(out,na.rm=TRUE))
}
