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
    out <- do.call(cbind.DNAbin, x@dna[genes])
    return(out[x@labels,,drop=FALSE])
}
