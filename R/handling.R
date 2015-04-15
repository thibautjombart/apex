#'
#' Concatenate genes into a single matrix
#'
#' This function concatenates separate DNA alignments into a single DNAbin matrix.
#'
#' @param x a \linkS4class{multidna} object.
#' @param genes an optional vector indicating the genes to retain for the concatenation; any way to subset the list in x@@dna is acceptable; by default, all genes are used.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @aliases concatenate
#'
#' @export
#'
#' @examples
#'
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' plot(x)
#'
#' image(concatenate(x))
concatenate <- function(x, genes=TRUE){
    if(inherits(x, "multidna")){
        out <- do.call(cbind.DNAbin, x@dna[genes])
        return(out[x@labels,,drop=FALSE])
    }
    if(inherits(x, "multiphyDat")){
        out <- do.call(cbind.phyDat, x@dna[genes])
        return(out)
    }
    return(x)
}





#################
## '[' operator
#################
#' Subset multidna objects
#'
#' Individuals in a \linkS4class{multidna} object can be subsetted like the rows of a matrix, with the form x[i,].
#' Genes can be subsetted like the columns of a matrix, i.e. with the form x[,j].
#'
#' @aliases [,multidna-method
#' @aliases [.multidna
#'
#' @param x the \linkS4class{multidna} object to subset.
#' @param i a vector of logical, integers or characters to subset data by individuals; characters will be matched against individual labels.
#' @param j a vector of logical, integers or characters to subset data by genes; characters will be matched against gene names labels.
#' @param drop present for compatibility with the generic; currently not used.
#' @param ... further arguments to be passed to other methods; currently ignored.
#'
#' @docType methods
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @export
#'
#' @examples
#'
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#' plot(x)
#'
#' ## keep only the first 5 individuals
#' x[1:5,]
#' plot(x[1:5,])
#'
#' ## keep individuals 2,4,6 and the second gene
#' x[c(2,4,6),2]
#' plot(x[c(2,4,6),2])
#'
setMethod("[", signature(x="multidna", i="ANY", j="ANY", drop="ANY"), function(x, i, j, ...) {
    if (missing(i)) i <- TRUE
    if (missing(j)) j <- TRUE

    ## subset data
    for(k in 1:length(x@dna)) x@dna[[k]] <- x@dna[[k]][i,,drop=FALSE]
    x@labels <- x@labels[i]
    x@dna <- x@dna[j]

    ## adjust counters
    x@n.ind <- length(x@labels)
    x@n.seq <- sum(sapply(x@dna, nrow))
    x@n.seq.miss <- .nMissingSequences(x@dna)
    return(x)
}) # end [] for SNPbin
