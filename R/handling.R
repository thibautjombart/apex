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
concatenate <- function(x, genes=TRUE){
    out <- do.call(cbind.DNAbin, x@dna[genes])
    return(out)
}
