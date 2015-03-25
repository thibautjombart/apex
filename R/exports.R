#'
#' Convert multidna into genind
#'
#' This function concatenates separate DNA alignments, and then extracts SNPs of the resulting alignment into a \linkS4class{genind} object.
#'
#' @param x a \linkS4class{multidna} object.
#' @param genes an optional vector indicating the genes to retain for the concatenation; any way to subset the list in x@@dna is acceptable; by default, all genes are used.
#' @param mlst if \code{TRUE}, each gene will result in a single locus in the genind object. (Default to \code{FALSE})
#' @param gapIsNA if \code{TRUE} and \code{mlst = TRUE}, sequences that consist entirely of gaps will be considered as NAs. (Default to \code{FALSE})
#'
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}, Zhian N. Kamvar
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
#' @importFrom adegenet DNAbin2genind df2genind
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
multidna2genind <- function(x, genes=TRUE, mlst=FALSE, gapIsNA=FALSE){
  if (!mlst){
    return(DNAbin2genind(concatenate(x, genes=genes)))
  } 
  xlist  <- lapply(x@dna, function(i) apply(as.character(i), 1, paste, collapse = ""))
  xdf    <- data.frame(xlist)
  xlevs  <- lapply(xdf, levels)
  if (gapIsNA){
    xgap <- find_gap_sequence(xlevs)
    for (i in names(xgap)){
      levels(xdf[[i]])[xgap[[i]]] <- NA
    }    
  }
  xdfnum <- data.frame(lapply(xdf, as.numeric))
  xgid   <- df2genind(xdfnum, ploidy = 1, ind.names = x@labels)
  names(xlevs)   <- names(xgid@all.names)
  xgid@all.names <- xlevs
  xgid@other$ind.info <- x@ind.info
  xgid@other$gene.info <- x@gene.info
  return(xgid)
}

find_gap_sequence <- function(x){
  wheregaps <- lapply(x, function(i) which(grepl("^\\-+?$", i)))
  return(wheregaps)
}

