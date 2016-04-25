#' Generics for DBAmmd-Class
#'
#' Generics for DBAmmd-Class
#'
#' @include DBAmmd-Class.R
#'
#' @inheritParams DBAmmd-Accessors
#' @param ... additional parameters
#' @rdname myGeneric

#' @rdname myGeneric
#' @export
setGeneric("metaData", function(x,...) standardGeneric("metaData"))

#' @rdname myGeneric
#' @export
setGeneric("Regions", function(x,...) standardGeneric("Regions"))

#' @rdname myGeneric
#' @export
setGeneric("Reads", function(x,...) standardGeneric("Reads"))

#' @rdname myGeneric
#' @export
setGeneric("Counts", function(x,...) standardGeneric("Counts"))

#' @rdname myGeneric
#' @export
setGeneric("Hists", function(x,...) standardGeneric("Hists"))

#' @rdname myGeneric
#' @export
setGeneric("Dists", function(x,...) standardGeneric("Dists"))

#' @rdname myGeneric
#' @export
setGeneric("Contrast", function(x,...) standardGeneric("Contrast"))

###

#' @rdname myGeneric
#' @export
setGeneric("numPeaks", function(x,...) standardGeneric("numPeaks"))

#' @rdname myGeneric
#' @export
setGeneric("numSamples", function(x,...) standardGeneric("numSamples"))

#' @rdname myGeneric
#' @export
setGeneric("Samples", function(x,...) standardGeneric("Samples"))

#' @rdname myGeneric
#' @export
setGeneric("Genome", function(x,...) standardGeneric("Genome"))


####--------
## Defining the slot setters
####--------

#' @rdname myGeneric
#' @export
setGeneric("setRegions", function(x,...) standardGeneric("setRegions"))

#' @rdname myGeneric
#' @export
setGeneric("setContrast", function(x,...) standardGeneric("setContrast"))

