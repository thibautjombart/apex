
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
    if(is.character(i)) i <- as.integer(na.omit(match(i, x@labels)))

    x@labels <- x@labels[i]
    x@dna <- x@dna[j]
    for(k in 1:length(x@dna)){
        toKeep <- x@labels[x@labels %in% rownames(x@dna[[k]])]
        x@dna[[k]] <- x@dna[[k]][toKeep,,drop=FALSE]
    }

    ## get rid of empty genes
    x@dna <- x@dna[sapply(x@dna, nrow)>0]

    ## adjust counters
    x@n.ind <- length(x@labels)
    x@n.seq <- sum(sapply(x@dna, nrow))
    x@n.seq.miss <- .nMissingSequences(x@dna)
    return(x)
}) # end [] for multidna