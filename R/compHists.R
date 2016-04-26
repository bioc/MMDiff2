#' Compute Peak histograms
#'
#' This function computes histograms at pre-defined regions (peaks)
#' from mapped fragments, i.e. fragment counts at genomic position. Note,
#' in contrast to genomic coverage or density maps, this function uses a single
#' position per fragment (usually its center) rather than the whole extend of
#' the fragment. This results in a significant increase in resolution.
#' The parameter \code{whichPos} determines whether
#' fragment centers, start or end positions should be considered
#' ('Center','Left','Right').
#' Results are stored as a list in the \code{Hists} slot of the DBAmmd Object,
#' with one entry per peak. For each peak i, a (n x L_i) matrix is generated,
#' where n is the number of samples and L_i is the number of bins used to cover
#' the extend of the peak. Note, L_i varies between peaks of different lengths.
#'
#' @inheritParams getPeakReads
#' @inheritParams DBAmmd-Accessors
#' @param bin.length size of binning window (in bp)
#' (DEFAULT: 20)
#'
#' @return DBAmmd object with updated slot \code{Hists}
#'
#' @seealso \code{\link{DBAmmd}}, \code{\link{getPeakReads}},
#' \code{\link{estimateFragmentCenters}}, \code{\link{plotPeak}},
#'
#'
#' @examples
#'
#' ## Example using a small data set provided with this package:
#'
#' data("MMD")
#' bin.length <- 20
#' MMD.1 <- compHists(MMD,bin.length)
#'
#' # use \code{plotPeak()} to plot indivdual peaks:
#' Peak.id <- '6'
#' plotPeak(MMD.1, Peak.id=Peak.id)
#'
#' # or explicitly using the histograms:
#' H <- Hists(MMD.1, whichPos='Center')
#' Sample <- 'WT.AB2'
#' Peak.idx <- match(Peak.id, names(Regions(MMD.1)))
#' plot(H[[Peak.idx]][Sample,],t='l')
#'
#' # add peak cooridnates:
#' Peak <- Regions(MMD.1)[Peak.idx]
#' meta <- metaData(MMD.1)
#' PeakBoundary <- meta$AnaData$PeakBoundary
#' x.coords <- as.integer(colnames(H[[Peak.idx]])) + start(Peak) - PeakBoundary
#' plot(x.coords,H[[Peak.idx]]['WT.AB2',],t='l',
#'     xlab=names(H)[Peak.idx], ylab='counts', main=Sample)
#'
#'
#' @export

# subfunctions:
#         - compHists
#         - compPeakMatrix
#         - addHists
#
# Gabriele Schweikert
# January 2014


compHists <- function(MD, bin.length=20, whichPos="Center"){

  message('checking parameters...')

  if (missing(MD))
    stop("MD object")

  Meta <- metaData(MD)
  PeakBoundary <- Meta$AnaData$PeakBoundary
  Peaks <- Regions(MD)

  VALs=c('Left.p','Right.n','Left.n','Right.p',
         'Center.n','Center.p',
         'Center','Left','Right')
  if (class(bin.length)!='numeric')
    stop("bin.length must be numeric")
  if (any(!is.element(whichPos,VALs))){
    message('whichPos have to be elements of:\n')
    print(VALs)
    stop("wrong whichPos")
  }

  ###determine which histogramms to compute
  todo <- intersect(whichPos,c("Center.p","Left.p","Right.p",
                               "Center.n","Left.n","Right.n"))
  if (any(is.element(whichPos,"Center"))){
    todo=c(todo,"Center.p","Center.n")
  }
  if (any(is.element(whichPos,"Left"))){
    todo=c(todo,"Left.p","Left.n")
  }
  if (any(is.element(whichPos,"Right"))){
    todo=c(todo,"Right.p","Right.n")
  }

  todo <- unique(todo)
  if (Meta$AnaData$paired==FALSE){
    todo <- intersect(todo, c("Left.p","Right.n","Center.p","Center.n"))
  }

  PeakHists <-  vector("list", length(todo))
  names(PeakHists) <- todo
  PeakLengths <- width(Peaks) + PeakBoundary*2
  for (j in 1:length(todo)){
    message('\n starting with ', todo[j])
    Reads <- Reads(MD,todo[j])

    Hists <- vector("list", length(Reads))
    names(Hists) <- names(Reads)

    for (i in 1:length(Reads)){
      message('sample ', names(Reads)[i])
      reads <- Reads[[i]]
      Hists[[i]] <- compHist(reads,PeakLengths, bin.length = bin.length)
    }
    PeakHists[[j]] <- compPeakMatrix(Hists,'Counts')
  }
  if (is.element("Center",whichPos)){
    PeakHists$Center <- addHists(PeakHists$Center.p,PeakHists$Center.n)
  }
  if (is.element("Left",whichPos)){
    if  (Meta$AnaData$pairedEnd){
      PeakHists$Left <- addHists(PeakHists$Left.p,PeakHists$Left.n)
    } else {
      PeakHists$Left <- PeakHists$Left.p
    }
  }
  if (is.element("Right",whichPos)){
    if  (Meta$AnaData$pairedEnd){
      PeakHists$Right <- addHists(PeakHists$Right.p,PeakHists$Right.n)
    } else {
      PeakHists$Right <- PeakHists$Right.n
    }
  }
  ## prepare output

  Meta$AnaData$Hist.bin.length=bin.length
  MD@MetaData <- Meta
  MD@Hists <- PeakHists

  return(MD)
}



# compHist takes a vector of positions and computes
# histogramms.
#
# INPUT   - Pos: list of lengths N with reads
#                  mapping to N peaks
#
#         - bin.length = 20 :
#                  number of base pairs that are collated into one bin
#
#         - Boundary = 200 : This is only
#                  used to determine start and end position of the
#                  histogram as 5'end of reads can be outside of peak
#                  region, histogram therefore covers the region: peak
#                  start - Boundary ...peak end + Boundary
#                  shift
#
# OUTPUT   - Hists$Counts: list of histograms of counts (length=N)
#          - Hists$Mids: Mid points for histogramms
#
# Gabriele Schweikert
# January 2016

