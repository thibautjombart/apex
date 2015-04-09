
######################
####  CONSTRUCTOR ####
######################
#' multiphyDat constructor
#'
#' New \linkS4class{multiphyDat} objects can be created using \code{new("multiphyDat", ...)} where "..." are arguments documented below.
#' The main input is a list of phyDat matrices. The constructor ensures that all matrices will be reordered in the same way, and genes with missing individuals will be filled by sequences of gaps ("-").
#'
#' @author Klaus Schliep \email{klaus.schliep@@gmail.com}
#' @author Thibaut Jombart \email{t.jombart@@imperial.ac.uk}
#'
#' @export
#'
#' @aliases initialize,multiphyDat-methods new.multiphyDat
#'
#' @param .Object the object skeleton, automatically generated when calling \code{new}.
#' @param dna a list of phyDat matrices (1 per gene); rows should be labelled and indicate individuals, but different individuals and different orders can be used in different matrices.
#' @param ind.info an optional data.frame containing information on the individuals, where individuals are in rows.
#' @param gene.info an optional data.frame containing information on the genes, where genes are in rows.
#' @param quiet a logical indicating if messages should be shown; defaults to FALSE.
#' @param ... further arguments to be passed to other methods
#'
#' @seealso
#' \itemize{
#' \item the \linkS4class{multiphyDat} class
#' \item \code{\link{read.multiphyDat}}
#' }
#' @examples
#' data(Laurasiatherian)
#' #' ## empty object
#' new("multiphyDat")
#'
#' ## simple conversion with nicely ordered output
#' \dontrun{
#' genes <- list(gene1=subset(Laurasiatherian,, 1:1600, FALSE),
#'     gene2=subset(Laurasiatherian,, 1601:3179, FALSE))
#' x <- new("multiphyDat", genes)
#' x
#' }
#'
#' ## trickier conversion with missing sequences / wrong order
#' genes <- list(gene1=subset(Laurasiatherian, 1:40),
#'     gene2=subset(Laurasiatherian, 8:47))
#' x <- new("multiphyDat", genes)
#' x
#'
setMethod("initialize", "multiphyDat", function(.Object, dna=NULL, ind.info=NULL, gene.info=NULL, quiet=FALSE, ...) {

    ## RETRIEVE PROTOTYPED OBJECT ##
    x <- .Object


    ## ESCAPE IF NO DATA ##
    if(is.null(dna)) return(x)


    ## HANDLE DNA ##
    ## cases where an multiphyDat is provided ##
    if(inherits(dna, "multiphyDat")){
        ind.info <- dna@ind.info
        gene.info <- dna@gene.info
        dna <- dna@dna
        dna@dna <- NULL
        invisible(gc())
    }

    ## cases where no info provided ##
    if(is.null(dna)) return(x)
    if(is.matrix(dna)) dna <- list(dna)

    ## coerce items in DNA to matrices ##
    # dna <- lapply(dna, as.matrix)
    fun <- function(x)ifelse(is.matrix(x),nrow(x),length(x))
    N.SEQ <- sum(sapply(dna, fun))
    if(N.SEQ==0){
        x@dna <- NULL
        x@ind.info <- x@gene.info <- NULL
        return(x)
    }

    ## convert matrices of characters into phyDat ##
    N.GENES <- length(dna)  
    for(i in 1:N.GENES){
        if(is.character(dna[[i]])) dna[[i]] <- phyDat(dna[[i]])
    }

    ## replace with generic names if needed ##
    if(is.null(names(dna))) names(dna) <- paste("gene", 1:N.GENES, sep=".")


    ## AUXILIARY FUNCTIONS ##

    ## HANDLE LABELS ##
    ## handle missing labels ##
#    missing.labels <- any(sapply(dna, function(e) is.null(labels(e))))
#    if(missing.labels){
#        if(!quiet) message("[multiphyDat constructor] missing/incomplete labels provided - using generic labels.\n")
        ## error if varying numbers of rows
#        if(length(unique(sapply(a, nrow)))>1) stop("[multiphyDat constructor] no labels provided and varying number of sequences across genes - cannot assume individuals are identical.")
#        labels <- paste("individual", 1:nrow(dna[[1]]), sep=".")
#        for(i in 1:N.GENES) rownames(dna[[i]]) <- labels
#    }

    ## get list of all labels ##
    all.labels <- unique(unlist(lapply(dna, names)))
    N.IND <- length(all.labels)


    ## COMPLETE/SORT MATRICES OF DNA ##
#    dna <- lapply(dna, form.dna.matrix, all.labels)


    ## PROCESS META INFO ##
    ## ind.info
    if(!is.null(ind.info)){
        if(nrow(ind.info)>N.IND && !quiet) warning("[multiphyDat constructor] ind.info has more rows than there are individuals")
        ind.info <- as.data.frame(ind.info)
    }

    ## gene.info
    if(!is.null(gene.info)){
        if(nrow(gene.info)>N.GENES && !quiet) warning("[multiphyDat constructor] gene.info has more rows than there are genes")
        gene.info <- as.data.frame(gene.info)
    }


    ## FORM FINAL OUTPUT ##
    x@dna <- dna
    x@labels <- all.labels
    x@n.ind <- N.IND
    x@n.seq <- N.SEQ
    x@ind.info <- ind.info
    x@gene.info <- gene.info

    return(x)
}) # end multiphyDat constructor
