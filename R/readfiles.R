
#' Read multiple DNA alignments
#'
#' These functions read multiple DNA alignments and store the output in a \linkS4class{multidna} object.
#' They are relying on ape's original functions \code{\link[ape]{read.dna}} and \code{\link[ape]{read.FASTA}}.
#'
#' @rdname readfiles
#' @aliases read.multidna
#' @aliases read.multiFASTA
#'
#' @param files a vector of characters indicating the paths to the files to read from
#' @param ... further arguments passed to \code{\link[ape]{read.dna}}
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @seealso  \code{\link[ape]{read.dna}}, \code{\link[ape]{read.FASTA}}
#'
#'
#'
read.multidna <- function(files, ...){
    out <- new("multidna", dna=lapply(files, read.dna, ...))
    return(out)
}


#'
#' @rdname readfiles
#'
read.multiFATSA <- function(files){
    out <- new("multidna", dna=lapply(files, read.FASTA))
    return(out)
}
