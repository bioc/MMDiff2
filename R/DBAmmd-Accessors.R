#' Extract data from DBAmmd objects
#'
#' This help file describes different ways to access the slots and values
#' contained in \code{\link{DBAmmd}} objects.
#'
#' @param x a DBAmmd Object. This Object can be created using \code{DBAmmd()}.
#' As MetaData it needs to contain the path to the data directory and the
#' name of a sampleSheet csv file. Additionally, the Regions of interests must
#' be provided as GRanges Object.
#' @param whichPos specifies which relative positions of mapped fragments
#' should to be considered.
#' Can be one of: 'Left.p', 'Right.p': Start and end positions of fragments
#' mapping to positive strand, respectively.
#' ('Right.p' is not available for single-end reads)
#' 'Left.n','Right.n': Start and end positions of fragments
#' mapping to negative strand, respectively.
#' ('Left.n' is not available for single-end reads)
#' Additionally inferred positions: 'Center.n','Center.p','Center','Left','Right'.
#' @param whichCounts can be 'T': total counts, or
#' 'p','n': counts of reads mapping to positive, negative strand, respectively
#' @param dist.method which method should be used to determine distances
#' between samples. Currently only Maximum Mean Discrepancy (MMD)
#' and Kolmogorov-Smirnov (KS) implemented.
#' (DEFAULT: 'MMD')
#' @param whichContrast (DEFAULT: 1)
#' @param Regions GRanges Object
#' @param contrast how to set contrast, either 'byCondition', or 'byTissue'.
#' contrast can also be manually set (see vignette for details.)
#'
#' @examples
#' data("MMD")
#'
#' Samples(MMD)
#' Genome(MMD)
#' numPeaks(MMD)
#' numSamples(MMD)
#' metaData(MMD)
#' R <- Regions(MMD)
#' Pos <- Reads(MMD, whichPos='Center')
#' C <- Counts(MMD)
#' H <- Hists(MMD, whichPos='Center')
#' D <- Dists(MMD)
#' C1 <- Contrast(MMD)
#'
#' @name DBAmmd-Accessors
#' @rdname DBAmmd-Accessors
#' @include DBAmmd-Class.R AllGenerics.R
NULL

#' @return \code{Genome(x)} returns
#' @rdname DBAmmd-Accessors
setMethod("Genome", "DBAmmd", function(x) x@MetaData$ExpData$genome)

#' @return \code{Samples(x)} returns
#' @rdname DBAmmd-Accessors
setMethod("Samples", "DBAmmd",function(x) x@MetaData$ExpData$samples)

#' @return \code{numPeaks(x)} returns
#' @rdname DBAmmd-Accessors
setMethod("numPeaks", "DBAmmd",function(x) length(x@rowRanges))

#' @return \code{numSamples(x)} returns
#' @rdname DBAmmd-Accessors
setMethod("numSamples", "DBAmmd",function(x) dim(Samples(x))[1])

#' @return \code{metaData(x)} returns
#' @rdname DBAmmd-Accessors
setMethod("metaData", "DBAmmd",function(x) x@MetaData)

#' @return \code{Regions(x)} returns
#' @rdname DBAmmd-Accessors
setMethod("Regions", "DBAmmd", function(x) x@rowRanges)

#' @return \code{Reads(x,whichPos)} returns
#' @rdname DBAmmd-Accessors
setMethod("Reads", "DBAmmd", function(x,whichPos=NULL){
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

#' @return \code{Counts(x,whichCounts)} returns
#' @rdname DBAmmd-Accessors
setMethod("Counts", "DBAmmd", function(x,whichCounts='T') {
  if (whichCounts=='T') x@RawTotalCounts
  else if (whichCounts=='p') x@RawCounts.p
  else if (whichCounts=='n') x@RawCounts.n
  else
    stop("wrong whichCounts; should be one of 'T': Total,
         'p': positive strand,'n': negative strand")
})

#' @return \code{Hists(x,whichPos)} returns
#' @rdname DBAmmd-Accessors
setMethod("Hists", "DBAmmd", function(x,whichPos){
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

#' @return \code{Dists(x,dist.method)} returns
#' @rdname DBAmmd-Accessors
setMethod("Dists", "DBAmmd", function(x,dist.method=NULL) {
  D <- x@DISTs
  if (is.null(dist.method)) return(D)
  else d <- D[[dist.method]]
  return(d)
})

#' @return \code{Contrast(x,whichContrast)} returns
#' @rdname DBAmmd-Accessors
setMethod("Contrast", "DBAmmd", function(x,whichContrast=1)
  x@Contrasts[whichContrast])


###SLOT SETTERs

#' @return \code{setRegions(x,Regions)} returns a DBAmmd Object with Regions set
#' @rdname DBAmmd-Accessors
setMethod("setRegions", "DBAmmd", function(x,Regions) {
  if (length(x@rowRanges)>0)
    warning('Object already contains rowRegions, will be overwritten.
            Check for inconsistencies')
  x@rowRanges=Regions;
  return(x)})

#' @return \code{setContrast(x,contrast)} returns a DBAmmd Object
#' with a contrast set
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



