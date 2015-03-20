#'
#' Convert multidna into genind
#'
#' This function concatenates separate DNA alignments, and then extracts SNPs of the resulting alignment into a \linkS4class{genind} object.
#'
#' @param x a \linkS4class{multidna} object.
#' @param genes an optional vector indicating the genes to retain for the concatenation; any way to subset the list in x@@dna is acceptable; by default, all genes are used.
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @aliases multidna2genind
#'
#' @seealso
#' \itemize{
#' \item concatenate
#' \item \code{\link{DNAbin2genind}} to convert single DNAbin objects.
#' }
#'
#' @export
#'
#' @importFrom adegenet DNAbin2genind
#'
#' @return a \linkS4class{genind} object
#'
#' @examples
#'
#' ## simple conversion with nicely ordered output
#' data(woodmouse)
#' genes <- list(gene1=woodmouse[,1:500], gene2=woodmouse[,501:965])
#' x <- new("multidna", genes)
#' x
#'
#'
multidna2genind <- function(x, genes=TRUE){
    return(DNAbin2genind(concatenate(x, genes=genes)))
}

