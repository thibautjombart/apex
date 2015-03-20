
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
#' @seealso
#' \itemize{
#' \item \code{\link[ape]{read.dna}}
#' \item  \code{\link[ape]{read.FASTA}}
#' }
#'
#' @export
#'
#' @examples
#' ## get path to the files
#' files <- dir(system.file(package="apex"),patter="patr", full=TRUE)
#' files
#'
#' ## read files
#' x <- read.multiFATSA(files)
#' x
#' plot(x)
#'
read.multidna <- function(files, ...){
    gene.names <- gsub(".fasta","",sapply(strsplit(files, "/"), tail, 1))
    dna <- lapply(files, read.dna, ...)
    names(dna) <- gene.names
    out <- new("multidna", dna=dna)
    return(out)
}


#'
#' @rdname readfiles
#' @export
read.multiFATSA <- function(files){
    gene.names <- gsub(".fasta","",sapply(strsplit(files, "/"), tail, 1))
    dna <- lapply(files, read.FASTA)
    names(dna) <- gene.names
    out <- new("multidna", dna=dna)
    return(out)
}
