#' Extract data from DBAmmd objects
#'
#' This help file describes different ways to access the slots and values
#' contained in a \code{\link{DBAmmd-class}} objects.
#'
#' @param x a DBAmmd Object. An empty instance can be created using \code{DBAmmd()}.
#' (See \code{\link{DBAmmd-class}} for more details.)
#' @param whichPos specifies which relative positions of mapped fragments
#' should to be considered.
#' Can be one of: 'Left.p', 'Right.p', 'Right.p' and 'Left.n':
#' Start and end positions of fragments mapping to positive or negative strand,
#' respectively ('Right.p' and 'Left.n' are not available for single-end reads).
#' Additionally inferred positions: 'Center.n','Center.p','Center','Left','Right'.
#' (DEFAULT: 'Center')
#' @param whichCounts can be 'T': total counts, or
#' 'p','n': counts of reads mapping to positive, negative strand, respectively.
#' @param dist.method specify method used for distances
#' between samples. Currently only Maximum Mean Discrepancy (MMD)
#' and Kolmogorov-Smirnov (KS) implemented.
#' (DEFAULT: 'MMD')
#' @param whichContrast index determining which of the set contrast should be used.
#' (DEFAULT: 1)
#' @param Regions GRanges Object specifying the Regions of Interesst / Peaks.
#' @param contrast determines how to set a new contrast for differential analysis.
#' A contrast can be automatically created either 'byCondition', or 'byTissue'.
#' The Contrast can also be manually set (see vignette for details).
#'
#' @seealso \code{\link{DBAmmd-class}}
#' @examples
#' data("MMD")
#'
#' Samples(MMD)
#' Genome(MMD)
#' numPeaks(MMD)
#' numSamples(MMD)
#' metaData(MMD)
#' R <- Regions(MMD)
#' Pos <- Reads(MMD)
#' C <- Counts(MMD)
#' H <- Hists(MMD)
#' D <- Dists(MMD)
#' C1 <- Contrast(MMD)
#'
#' @name DBAmmd-Accessors
#' @rdname DBAmmd-Accessors
#' @include DBAmmd-Class.R AllGenerics.R
NULL

#' @return \code{Genome(x)} returns the name of the used genome version, if set
#' in the metaData.
#' @rdname DBAmmd-Accessors
setMethod("Genome", "DBAmmd", function(x) x@MetaData$ExpData$genome)

#' @return \code{Samples(x)} returns the information which was provided in the
#' SampleSheet.csv to describe the data.
#' @rdname DBAmmd-Accessors
setMethod("Samples", "DBAmmd",function(x) x@MetaData$ExpData$samples)

#' @return \code{numPeaks(x)} returns the number of Peaks / Regions of Interest
#' that are associated with the DBAmmd object.
#' @rdname DBAmmd-Accessors
setMethod("numPeaks", "DBAmmd",function(x) length(x@rowRanges))

#' @return \code{numSamples(x)} returns the number of samples associated with the
#' DBAmmd object.
#' @rdname DBAmmd-Accessors
setMethod("numSamples", "DBAmmd",function(x) dim(Samples(x))[1])

#' @return \code{metaData(x)} returns the metaData associated with the
#' DBAmmd object.
#' @rdname DBAmmd-Accessors
setMethod("metaData", "DBAmmd",function(x) x@MetaData)

#' @return \code{Regions(x)} returns the Peaks / Regions of Interest that are
#' associated with the DBAmmd object.
#' @rdname DBAmmd-Accessors
setMethod("Regions", "DBAmmd", function(x) x@rowRanges)

#' @return \code{Reads(x,whichPos)} returns the Reads mapping to the Regions of Interest.
#' @rdname DBAmmd-Accessors
setMethod("Reads", "DBAmmd", function(x,whichPos='Center'){
  reads=x@Reads

  VALs=c('Left.p','Right.n','Left.n','Right.p','Center.n','Center.p',
         'Center','Left','Right')
  if (is.null(whichPos) || !is.element(whichPos,VALs)){
    message('whichPos has to be one of:\n')
    print(VALs)
    stop("wrong whichPos'")
  }
  R=reads[[whichPos]]
  return(R)
})

#' @return \code{Counts(x,whichCounts)} returns a m x n matrix containing the
#' Counts of Reads mapping to the Peaks / Regions of Interest.
#' Depending on the value of 'whichCounts', total counts ('T'),
#' or counts of reads mapping to positive ('p'), or negative strand ('n')
#' are returnt. See \code{\link{getPeakReads}} for more details.
#' @rdname DBAmmd-Accessors
setMethod("Counts", "DBAmmd", function(x,whichCounts='T') {
  if (whichCounts=='T') x@RawTotalCounts
  else if (whichCounts=='p') x@RawCounts.p
  else if (whichCounts=='n') x@RawCounts.n
  else
    stop("wrong whichCounts; should be one of 'T': Total,
         'p': positive strand,'n': negative strand")
})

#' @return \code{Hists(x,whichPos)} returns a list of matrices of length m
#' (number of Peaks). Each matrix is a n x L_i matrix, where n is the number of
#' samples and L_i is the number of bins used to cover
#' the extend of the peak. Note, L_i varies between peaks of different lengths.
#' See \code{\link{compHists}} for more details.
#' @rdname DBAmmd-Accessors
setMethod("Hists", "DBAmmd", function(x,whichPos='Center'){
  hists=x@Hists
  VALs=c('Left.p','Right.n','Left.n','Right.p','Center.n','Center.p',
         'Center','Left','Right')
  H=hists[[whichPos]]
  if (!is.element(whichPos,VALs)){
    message('whichPos has to be one of:\n')
    print(VALs)
    stop("wrong whichPos'")
  }
  return(H)
}
)

#' @return \code{Dists(x,dist.method)} returns a matrix containing distances
#' between pairs of samples for each peak. See \code{\link{compDists}} for
#' more details.
#' @rdname DBAmmd-Accessors
setMethod("Dists", "DBAmmd", function(x,dist.method=NULL) {
  D <- x@DISTs
  if (is.null(dist.method)) return(D)
  else d <- D[[dist.method]]
  return(d)
})

#' @return \code{Contrast(x,whichContrast)} returns the specified contrast.
#' @rdname DBAmmd-Accessors
setMethod("Contrast", "DBAmmd", function(x,whichContrast=1)
  x@Contrasts[whichContrast])


###SLOT SETTERs

#' @return \code{setRegions(x,Regions)} returns a DBAmmd Object with set
#' Peaks / Regions of Interests.
#' @rdname DBAmmd-Accessors
setMethod("setRegions", "DBAmmd", function(x,Regions) {
  if (length(x@rowRanges)>0)
    warning('Object already contains rowRegions, will be overwritten.
            Check for inconsistencies')
  x@rowRanges=Regions;
  return(x)})

#' @return \code{setContrast(x,contrast)} returns a DBAmmd Object
#' with a set contrast.
#' @rdname DBAmmd-Accessors
setMethod("setContrast", "DBAmmd", function(x,contrast) {
  if (contrast=='byCondition' | contrast=='byTissue'){
    if (contrast=='byCondition'){
      by.what = unique(Samples(x)$Condition)
      what = 'Condition'
    } else if (contrast=='byTissue'){
      by.what = unique(Samples(x)$Tissue)
      what = 'Tissue'
    }
    if (length(by.what)!=2){
      stop('need 2 conditions in Samplesheet to set contrast;
           set contrast by manually')
    }
    group1 <- Samples(x)[[what]]==by.what[1]
    group2 <- Samples(x)[[what]]==by.what[2]
    names(group1) <- Samples(x)$SampleID
    names(group2) <-  Samples(x)$SampleID
    contrast <- list(group1=group1,
                     group2=group2,
                     name1=by.what[1],
                     name2=by.what[2])
    } else if (class(contrast)!=list){
      stop('wrong contrast')
    }

  contrasts <- x@Contrasts
  if (length(contrasts)>0){
    message('Object already contains contrasts,
            current contrast will be appended')
    contrasts[[length(contrasts)+1]] <- contrast
  } else{
    contrasts <- list(contrast)
  }
  x@Contrasts <- contrasts
  return(x)})



