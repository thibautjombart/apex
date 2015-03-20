







## ####################
## ####  ACCESSORS ####
## ####################

## ################
## ## get.nlocus ##
## ################
## setMethod("get.nlocus","multidna", function(x, ...){
##     if(is.null(x@dna)) return(0)
##     return(length(x@dna))
## })



## ################
## ## get.locus ##
## ################
## setMethod("get.locus","multidna", function(x, ...){
##     if(is.null(x@dna)) return(NULL)
##     return(names(x@dna))
## })



## ###################
## ## get.sequences ##
## ###################
## ##  (get sequence IDs)
## setMethod("get.sequences","multidna", function(x, ...){
##     if(is.null(x)) return(NULL)
##     return(unlist(lapply(x@dna, rownames)))
## })



## ####################
## ## get.nsequences ##
## ####################
## setMethod("get.nsequences","multidna", function(x, what=c("total","bylocus"), ...){
##     what <- match.arg(what)
##     nLoc <- get.nlocus(x)
##     if(nLoc==0) return(0)

##     temp <- sapply(x@dna, nrow)
##     if(what=="bylocus") return(temp)
##     return(sum(temp))
## })



## #####################
## ## get.individuals ##
## #####################
## setMethod("get.individuals","multidna", function(x, ...){
##     if(is.null(x)) return(NULL)
##     return(unique(x@meta$individualID))
## })



## ######################
## ## get.nindividuals ##
## ######################
## setMethod("get.nindividuals","multidna", function(x, ...){
##     if(is.null(x)) return(0)
##     return(length(get.individuals(x)))
## })



## ###############
## ## get.dates ##
## ###############
## setMethod("get.dates","multidna", function(x, ...){
##     if(is.null(x)) return(NULL)
##     return(unique(x@meta$date))
## })



## ################
## ## get.ndates ##
## ################
## setMethod("get.ndates","multidna", function(x, ...){
##     if(is.null(x)) return(0)
##     return(length(get.dates(x)))
## })



## #############
## ## get.dna ##
## #############
## ## returns a matrix of dna sequences for a given locus
## setMethod("get.dna","multidna", function(x, locus=NULL, id=NULL, ...){
##     ## return NULL if no info ##
##     nLoc <- get.nlocus(x)
##     if(nLoc==0) return(NULL)

##     ## RETURN SLOT CONTENT AS IS IF NOTHING ELSE ASKED ##
##     if(is.null(locus) && is.null(id)) return(x@dna)

##     ## INFO REQUESTED PER LOCUS ##
##     if(is.null(id)){
##         ## return only locus if nLoc==1 and no info on locus ##
##         if(nLoc==1 && is.null(locus)) return(x@dna[[1]])

##         ## otherwise use locus info ##
##         if(nLoc>1 && is.null(locus)) stop("locus must be specified (data contain more than one locus)")
##         return(x@dna[locus])
##     }

##     ## INFO REQUESTED PER SEQUENCE ID ##
##     ## if logicals or integers, find corresponding names
##     if(is.logical(id) | is.numeric(id) | is.integer(id)){
##         id <- get.sequences(x)[id]
##     }
##     id <- as.character(id)
##     if(!all(id[!is.na(id)] %in% get.sequences(x))) {
##         temp <- paste(id[!is.na(id) & !id %in% get.sequences(x)], collapse=", ")
##         warning(paste("The following sequence IDs are not in the dataset:", temp))
##         id <- id[!is.na(id) & id %in% get.sequences(x)]
##     }
##     out <- lapply(x@dna, function(e) e[id[id %in% rownames(e)],,drop=FALSE])
##     out <- out[sapply(out, nrow)>0]
##     return(out)
## })







