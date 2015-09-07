## ####################
## ####  ACCESSORS ####
## ####################
#' @name accessors
#' @title multidna Accessors
#' @description Accessors for slots in \linkS4class{multidna} objects.
#'
#' @param x a \linkS4class{multidna} object.
#' @param loci a character, numeric, or logical vector identifying which
#'   loci to return.
#' @param ids a character, numeric, or logical vector identifying which
#'   sequences to return within a locus.
#' @param gap.only logical. Return information only for sequences containing
#'   all gaps?
#' @param simplify logical. If \code{FALSE}, always return a list of
#'   DNAbin sequences. If \code{TRUE} and only one locus has been requested,
#'   return a single DNAbin object.
#' @param exclude.gap.only logical. Remove sequences containing all gaps?
#' @param value a replacement value for the slot.
#' @param ... further arguments passed on to other functions.
#'
#' @details
#' \describe{
#'   \item{getNumLoci}{Returns the number of loci.}
#'   \item{locusNames}{Returns or sets the names of each locus.}
#'   \item{getNumSequences}{Returns the number of sequences in each locus.}
#'   \item{getSequenceNames}{Returns the names of individual sequences at each
#'     locus.}
#'   \item{getSequences}{Returns sequences of specified loci and individuals.}
#' }
#'
#setClass("multidna")

## ################
## ## getNumLoci ##
## ################
#' @rdname accessors
#' @export
setGeneric("getNumLoci", function(x, ...) standardGeneric("getNumLoci"))
#' @rdname accessors
#'
#' @aliases getNumLoci
#' @aliases getNumLoci,multidna
#' @export
setMethod("getNumLoci", "multidna", function(x, ...) {
  if(is.null(x@dna)) return(0)
  return(length(x@dna))
})


## ################
## ## locusNames ##
## ################
#' @rdname accessors
#' @export
setGeneric("locusNames", function(x, ...) standardGeneric("locusNames"))
#' @rdname accessors
#' @aliases locusNames,multidna
#' @export
setMethod("locusNames", "multidna", function(x, ...) names(x@dna))
#' @rdname accessors
#' @export
setGeneric("locusNames<-", function(x, value) standardGeneric("locusNames<-"))
#' @rdname accessors
#' @aliases locusNames<-,multidna
#' @export
setMethod("locusNames<-", "multidna", function(x, value) {
  names(x@dna) <- value
  validObject(x)
  x
})


## #####################
## ## getNumSequences ##
## #####################
#' @rdname accessors
#' @export
setGeneric("getNumSequences", function(x, ...) standardGeneric("getNumSequences"))
#' @rdname accessors
#' @aliases getNumSequences,multidna
#' @export
setMethod("getNumSequences", "multidna",
          function(x, gap.only = FALSE, loci = NULL, ...) {
  # check that object isn't empty
  if(is.null(x@dna)) {
    warning("'x' is empty. NULL returned.", call. = FALSE)
    return(NULL)
  }
  loci <- .checkLocusNames(x, loci)

  sapply(loci, function(this.locus) {
    dna <- x@dna[[this.locus]]
    if(gap.only) sum(.isGapOnly(dna)) else nrow(as.matrix(dna))
  })
})


## ######################
## ## getSequenceNames ##
## ######################
#' @rdname accessors
#' @export
setGeneric("getSequenceNames", function(x, ...) standardGeneric("getSequenceNames"))
#' @rdname accessors
#' @aliases getSequenceNames,multidna
#' @export
setMethod("getSequenceNames", "multidna",
          function(x, gap.only = FALSE, loci = NULL, ...) {
  # check that object isn't empty
  if(is.null(x@dna)) {
    warning("'x' is empty. NULL returned.", call. = FALSE)
    return(NULL)
  }
  loci <- .checkLocusNames(x, loci)

  sapply(loci, function(this.locus) {
    dna <- x@dna[[this.locus]]
    if(gap.only) labels(dna)[.isGapOnly(dna)] else labels(dna)
  }, simplify = FALSE)
})


## ##################
## ## getSequences ##
## ##################
#' @rdname accessors
#' @export
setGeneric("getSequences", function(x, ...) standardGeneric("getSequences"))
#' @rdname accessors
#' @aliases getSequences,multidna
#' @export
setMethod("getSequences", "multidna",
          function(x, loci = NULL, ids = NULL, simplify = TRUE,
                   exclude.gap.only = TRUE, ...) {
  # check that object isn't empty
  if(is.null(x@dna)) {
    warning("'x' is empty. NULL returned.", call. = FALSE)
    return(NULL)
  }

  # loop through loci
  loci <- .checkLocusNames(x, loci)
  new.dna <- sapply(loci, function(this.locus) {
    # extract this DNAbin object
    dna <- as.list(x@dna[[this.locus]])
    if(exclude.gap.only) dna <- dna[!.isGapOnly(dna)]
    # return sequences for IDs which are present
    locus.ids <- .checkIDs(dna, ids)
    if(is.null(locus.ids)) NULL else dna[locus.ids]
  }, simplify = FALSE)
  new.dna <- new.dna[!sapply(new.dna, is.null)]

  if(length(new.dna) == 1 & simplify) new.dna[[1]] else new.dna
})







