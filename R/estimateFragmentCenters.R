#' estimate center of fragments
#'
#' This function computes average shifts between
#' forward and reverse strand and applies it to estimate fragment centers.
#'
#' @inheritParams getPeakReads
#' @param shift can be set if the offset between forward and reverse strand is
#' known (e.g. 1/2 median fragment size). In this case shift will not be
#' estimated  (DEFAULT: NULL)
#' @param draw.on plot scatterplots for counts on
#' forward vs reverse strand and histograms of determined shifts
#' (DEFAULT: FALSE)
#'
#' @return DBAmmd object with updated slots \code{Reads} and \code{MetaData}.
#'
#' @examples
#'
#' ## Example using a small data set provided with this package
#'
#' data("MMD")
#' MMD.1 <- estimateFragmentCenters(MMD)
#'
#' # To access centers of fragments:
#' Reads.C <- Reads(MMD.1,'Center')
#'
#' # To access the determined shifts for each sample:
#' meta <- metaData(MMD.1)
#' meta$AnaData$Shifts
#'
#'
#' @seealso \code{\link{DBAmmd}}, \code{\link{getPeakReads}},
#' \code{\link{compHists}}
#'
#' @import graphics
#' @importFrom stats ccf median
#' @export

#
#
# subfunctions:
#         - getPeakIds4shift
#         - getShift
#         - shiftReads
#
#
# Gabriele Schweikert
# January 2014
# needs paired end!! checking cluster usage, speed-ups + memory improvement


estimateFragmentCenters <- function(MD, shift=NULL, draw.on=FALSE){

  message('checking parameters...')

  if (missing(MD))
    stop("MD object")

  Meta <- metaData(MD)
  PeakBoundary <- Meta$AnaData$PeakBoundary
  Peaks <- Regions(MD)

  Left.p <- Reads(MD,'Left.p')
  Right.n <- Reads(MD, 'Right.n')

  numSamples <-  length(Left.p)
  SampleIDs <- names(Left.p)
  Center.p <- vector("list", numSamples)
  names(Center.p) <- names(Left.p)
  Center.n <- Center.p
  Center <- Center.p

  if (!Meta$AnaData$pairedEnd){
    Shifts <- rep(0, numSamples)
    names(Shifts) <- SampleIDs
    for (i in 1:numSamples){
      ## CORRECT FOR SHIFT BETWEEN +/-STRAND

      Idx <- getPeakIds4shift(Counts(MD,'p')[,i],Counts(MD,'n')[,i], draw.on,
                              name = SampleIDs[i])
      if (is.null(shift)){
        if (length(Idx)==0){
          message('can\'t determine shift between forward and
                  reverse strand. Not enough peaks!')
        } else{
          PeakLengths <- width(Peaks)[Idx] +PeakBoundary*2
          Shifts[i] <- getShift(Left.p[[i]][Idx],Right.n[[i]][Idx],
                                PeakLengths, draw.on = draw.on,
                                name=SampleIDs[i])
        }
      } else  Shifts[i] <- shift

      message('estimateFragmentCenter')
      R <- list(Left.p=Left.p[[i]],Right.n = Right.n[[i]])
      R <- shiftReads(R, Shifts[i])

      Center.p[[i]] <- R$Center.p
      Center.n[[i]] <- R$Center.n
      Center[[i]] <- R$Center
      message("")
    } # end loop over files
  } else if (Meta$pairedEnd) {
    message('Paired end')
  }

  ## prepare output

  R <- list(Left.p = Left.p,
            Right.n = Right.n,
            Center.p = Center.p,
            Center.n = Center.n,
            Center = Center,
            Left = Left.p,
            Right = Right.n)


  Meta$AnaData$Shifts <- Shifts

  MD@MetaData <- Meta
  MD@Reads <- R
  return(MD)
}

######################
# getPeakIds4shifts: Determines which Peaks to use to for strand shift
# correction. Peaks are selected wich have total counts on positve and
# negative strand in the 9th 10-quantiles of all total counts.
#
# INPUT   - Counts.p  (nb of reads mapping to peak on + strand)
#         - Counts.n  (nb of reads mapping to peak on + strand)
#         - draw.on=TRUE
#         - name (for figure caption)
#
# OUTPUT  - Idx: index of peaks
#
#
# Gabriele Schweikert
# January 2014

getPeakIds4shift <- function(Counts.p,Counts.n,draw.on=TRUE,name=NULL){

  summary(Counts.p)
  summary(Counts.n)

  x <- quantile(Counts.p, probs=seq(0,1,0.1))
  ii1 <- which(Counts.p>x[9] & Counts.p<x[10])
  x <- quantile(Counts.n, probs=seq(0,1,0.1))
  ii2 <- which(Counts.n>x[9] & Counts.n<x[10])
  Idx <- intersect(ii2, ii1)

  message('Using ',length(Idx), ' peaks to determine shift')
  if (draw.on){
    if (is.null(name)){
      main = ''
    } else{
      main = name
    }
    dev.new()
    ## fn <- paste(pics.dir,'forwardVSReverse_',main,'pdf',sep='')
    ## pdf(fn)
    smoothScatter(Counts.p, Counts.n, xlab='total counts on + strand',
                  ylab='total counts on - strand', main=main)
    if (length(Idx)>0){
      points(Counts.p[Idx], Counts.n[Idx], col='red')
    }
    legend("topleft", pch=c(21),col='red', "peaks to determine strand shift")
    # dev.off()
  }
  return(Idx)
}



######################
# getshift: Determines strand shift between reads mapping to positive
# and negative strand
#
# INPUT   - Reads: Reads$Left.p
#                  Reads$Right.n
#         - Idx: index of peaks to use
#         - bin.length = 10
#         - draw.on = TRUE
#
# OUTPUT  - shift
#
#
# Gabriele Schweikert
# January 2014

getShift <- function(p,n, PeakLengths, bin.length=10, draw.on, name ){
  message('using bin.length = ', bin.length, ' to determine shift.')

  H.p <- compHist(p, PeakLengths, bin.length=bin.length)
  H.n <- compHist(n, PeakLengths, bin.length=bin.length)
  # H.p <- compHist(p, bin.length=bin.length)
  # H.n <- compHist(n, bin.length=bin.length)

  shift <- rep(0, length(H.p$Counts))
  for (i in 1:length(H.p$Counts)){
    a <- ccf(H.n$Counts[[i]], H.p$Counts[[i]], lag.max=40, plot=FALSE);
    if (length(which.max(a$acf))>0){
      shift[i] <- a$lag[which.max(a$acf)]
    } else {
      shift[i] <- 0
    }
  }

  shift <- shift*bin.length
  message('shifts statistics:')
  print(summary(shift))
  if (draw.on){
    if (is.null(name)){
      main = ''
    } else{
      main = name
    }
    dev.new()
    ## fn <- paste(picss.dir,'/strandShifts_',main,'pdf',sep='')
    ## pdf(fn)
    hist(shift, 100, xlab='shift between + and - strand [bp]', main=main)
    lines(x=c(median(shift), median(shift)), y=c(0,length(shift)),
          col='red', lwd=2)
    # dev.off()

  }
  shift <- median(shift)
  message('determined shift between positive and reverse strand is: ', shift )
  return(shift)
}

######################
# shiftReads: shifts Reads$p by +shift/2 and Reads$n by
# -shift/2
#
# INPUT   - Reads: Reads$p
#                  Reads$n
#         - shift
#
# OUTPUT  - Reads
#
#
# Gabriele Schweikert
# August 2012

shiftReads <- function(Reads, shift){

  p <- Reads$Left.p
  #n <- Reads$Left.n

  p <- lapply(p, function(P){P <- P+shift/2})
  #n <- lapply(n, function(P){P <- P-shift/2})

  Reads$Center.p <- p
  #Reads$Left.n <- n

  #p <- Reads$Right.p
  n <- Reads$Right.n

  #p <- lapply(p, function(P){P <- P+shift/2})
  n <- lapply(n, function(P){P <- P-shift/2})

  #Reads$Right.p <- p
  Reads$Center.n <- n

  Reads$Center <- mapply(c, p,n, SIMPLIFY=FALSE)

  return(Reads)
}


