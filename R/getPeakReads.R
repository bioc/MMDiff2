#' Get reads from indexed bam files for defined regions
#'
#' This function collects all short reads from bam files that map
#' to pre-defined regions of interest. Note, that it fetches the exact start and
#' end positions of mapped fragments, not the coverage. In the case of
#' single-end reads, the left most postions of fragments mapping to the positive
#' strands and the right most positions of fragments
#' mapping to the negative strands are stored. To find centers of fragments use
#' \code{estimateFragmentCenters()}. Positions are given relative to the start
#' of the peak. Also computed are TotalCounts, i.e. number of fragments mapping to a peak
#' region, as well as number of fragments mapping to forward and reverse strand.
#'
#' @param MD DBAmmd Object. This Object can be created using \code{DBAmmd()}.
#' @param PeakBoundary Defines flanking regions of peaks.
#' The estimated fragment length is a good guess
#' (DEFAULT: 200).
#' @param pairedEnd whether the reads are paired-end
#' (paired-end is currently not fully tested)
#' (DEFAULT: FALSE)
#' @param run.parallel whether to run in parallel
#' (currently no parallelization implemented)
#' (DEFAULT: FALSE)
#'
#'
#' @return DBAmmd object with updated slots
#'
#' @seealso \code{\link{DBAmmd}}, \code{\link{estimateFragmentCenters}}
#'
#' @examples
#'
#' ## Example using a small data set provided in the MMDiffBamSubset package
#'
#' # setting the Experiment meta data:
#' ExpData <- list(dataDir=system.file("extdata", package="MMDiffBamSubset"),
#'            sampleSheet="Cfp1.csv")
#'
#' MetaData <- list('ExpData' = ExpData)
#'
#' # Creating a DBAmmd data set:
#' MMD <- DBAmmd(MetaData)
#'
#' # defining a small Region for which to get reads:
#' Regions <- GRanges(seqnames=c('chr1'),
#'            IRanges(start = c(4560912,4677889), end = c(4562680,4679681)))
#' MMD <- setRegions(MMD,Regions)
#' MMD <- getPeakReads(MMD)
#'
#' # To access Left ends of fragments mapping to positive strand:
#' Reads.L <- Reads(MMD,'Left.p')
#'
#' # To access Right ends of fragments mapping to negative strand:
#' Reads.R <- Reads(MMD,'Right.n')
#'
#' # To access Matrix of TotalCounts:
#' C.t <- Counts(MMD,whichCounts='T')
#'
#' # Counts on positive strand:
#' C.p <- Counts(MMD,whichCounts='p')
#'
#' # Counts on negative strand:
#' C.n <- Counts(MMD,whichCounts='n')
#'
#' @import GenomicRanges Rsamtools grDevices parallel
#' @export
#
# subfunctions:
#         - getReadsWrapper
#         - getReads
#         - sortReads (needs to be up-dated for counts)
#
#
# Gabriele Schweikert
# March 2016
# needs checking for paired-end and cluster usage

getPeakReads <- function(MD,
                         PeakBoundary = 200,
                         pairedEnd = FALSE,
                         run.parallel = FALSE){



  message('checking parameters...')
  if (missing(MD))
    stop("DBAmmd object required")
  if (pairedEnd)
    message('assuming paired-end reads')
  if (!pairedEnd)
    message('assuming single-end reads')

  message('using Peak boundary:',PeakBoundary,"bp \n")

  if (!run.parallel)
    message('parallel running off')
  ## cluster stuff


  if (run.parallel){
    message('running parallel')
    NCPUs <- detectCores()
    cl <- makeCluster(NCPUs)
  }

  wd <- setwd(metaData(MD)$ExpData$dataDir)
  Peaks <- Regions(MD)
  chroms <- unique(as.character(seqnames(Peaks)))

  samples <- Samples(MD)
  numChips <- nrow(samples)
  chips <- as.matrix(samples$bamReads,ncol=1)
  rownames(chips) = samples$SampleID
  chips <- unique(chips)
  Todo <- chips
  if (!is.null(samples$bamControl)){
    inputs <- as.matrix(samples$bamControl,ncol=1)
    rownames(inputs) <- samples$ControlID
    inputs <- unique(inputs)
    if (length(inputs)>0){
      Todo <- rbind(chips, inputs)
    }
  }
  todo <- as.vector(Todo)
  names(todo) <- rownames(Todo)
  #check sample and control ids
  if (any(is.null(names(todo)))) {
    print(todo)
    stop('missing SampleIDs or ControlIDs')}
  if (anyDuplicated(names(todo))>0){
    print(todo)
    stop('SampleIDs, ControlIDs not unique')
  }
  ## message('Peak length statistics:')
  ## print(summary(width(Peaks)))
  FNS <- lapply(todo, function(FN){
    fn <- unlist(strsplit(FN,'/')); fn <- fn[length(fn)]})
  FNS <- lapply(FNS, function(FN){
    fn <- unlist(strsplit(FN,'.bam')); fn <- fn[1]})

  ## RawReads <- vector("list", length(todo))
  ## names(RawReads) <- FNS
  Left.p <- vector("list", length(todo))
  names(Left.p) <- names(todo)
  Right.n <- Left.p

  if (pairedEnd){
    Left.n <- Left.p
    Right.p <- Left.p
    Center.p <- Left.p
    Center.n <- Left.p
  }

  numPeaks <- length(Peaks)
  RawTotalCounts <- matrix(0,nrow=numPeaks,ncol=length(todo))
  colnames(RawTotalCounts) <- names(todo)
  RawCounts.p <- RawTotalCounts
  RawCounts.n <- RawTotalCounts

  for (i in 1:length(todo)){

    message('loading reads from file: ', todo[i])

    if (run.parallel){
      cl.idx <-  clusterSplit(cl,1:length(chroms))
      chrsplt <- lapply(cl.idx, function(id) chroms[id])
      r <- clusterApplyLB(cl,chrsplt, getReadsWrapper,Peaks,todo[[i]],
                          pairedEnd=pairedEnd)
      Reads <- sortReads(r,Peaks)
    } else {
      Reads <- getReadsWrapper(chroms,Peaks,todo[[i]], pairedEnd=pairedEnd)
    }

    ## Read positions are shifted to peak coordinates
    ## i.e. for each peak reads mapping exactly to the 5'end of the peak region
    ## are at pos = PeakBoundary+1

    starts <- start(Peaks)
    for (j in 1:length(starts)){
      Reads$Left.p[[j]] <- Reads$Left.p[[j]] - starts[j]+PeakBoundary+1
      Reads$Right.n[[j]] <- Reads$Right.n[[j]]  - starts[j]+PeakBoundary+1
      if (pairedEnd){
        Reads$Left.n[[j]] <- Reads$Left.n[[j]] - starts[j]+PeakBoundary+1
        Reads$Right.p[[j]] <- Reads$Right.p[[j]] - starts[j]+PeakBoundary+1
        Reads$Center.p[[j]] <- Reads$Center.p[[j]] - starts[j]+PeakBoundary+1
        Reads$Center.n[[j]] <- Reads$Center.n[[j]] - starts[j]+PeakBoundary+1
      }
    }
    Left.p[[i]] <- Reads$Left.p
    Right.n[[i]] <- Reads$Right.n
    if (pairedEnd){
      Left.n[[i]] <- Reads$Left.n
      Right.p[[i]] <- Reads$Right.p
      Center.p[[i]] <- Reads$Center.p
      Center.n[[i]] <- Reads$Center.n
    }

    RawTotalCounts[,i] <- Reads$RawTotalCounts
    RawCounts.p[,i] <- Reads$RawCounts.p
    RawCounts.n[,i] <- Reads$RawCounts.n
    message("")
  }# end loop over files

  ## prepare output

  if (pairedEnd){
    Reads <- list(Left.p = Left.p,
                  Left.n = Left.n,
                  Right.p = Right.p,
                  Right.n = Right.n,
                  Center.p = Center.p,
                  Center.n = Center.n)

  } else {
    Reads <- list(Left.p=Left.p,Right.n=Right.n)
  }

  Meta <- metaData(MD)
  Meta$AnaData <- list(pairedEnd=pairedEnd,
                       PeakBoundary=PeakBoundary)
  MD@MetaData <- Meta
  MD@Reads <- Reads
  MD@RawTotalCounts = RawTotalCounts
  MD@RawCounts.p=RawCounts.p
  MD@RawCounts.n=RawCounts.n


  if (run.parallel){
    stopCluster(cl)
  }
  setwd(wd)
  return(MD)
}


getReadsWrapper <- function(chroms,Peaks,bam.file,pairedEnd){

  Left.p <- vector("list", length(Peaks))
  Left.n <- Left.p
  Right.p <- Left.p
  Right.n <- Left.p
  Center.p <- Left.p
  Center.n <- Left.p
  RawTotalCounts <- rep(0,length(Peaks))
  RawCounts.p <- RawTotalCounts
  RawCounts.n <- RawTotalCounts
  peak.id <- 0
  if (length(chroms)==0){
    return()
  }

  for (j in 1:length(chroms)){
    if (is.element(chroms[j],seqnames(seqinfo(Peaks)))){
      P <- Peaks[seqnames(Peaks)==chroms[j]]
      message('-chromosome: ', chroms[j], ' containing ', length(P), ' peaks')
      P <- cbind(rep(j,length(P)),start(P),end(P))
      reads <- getReads(P, chroms[j],bam.file, pairedEnd=pairedEnd)
      ids <- peak.id+seq(1, nrow(P))

      if (length(ids) != length(reads$Left.p)){
        message('trouble on chromosome', chroms[j]) }

      Left.p[ids] <- reads$Left.p
      Right.n[ids] <- reads$Right.n
      names(Left.p)[ids] <- names(reads$Left.p)
      names(Right.n)[ids] <- names(reads$Right.n)
      RawTotalCounts[ids] <- reads$RawTotalCounts
      RawCounts.p[ids] <- reads$RawCounts.p
      RawCounts.n[ids] <- reads$RawCounts.n

      if (pairedEnd){
        Left.n <- c(Left.n,reads$Left.n)
        Right.p <- c(Right.p,reads$Right.p)
        Center.n <- c(Center.n,reads$Center.n)
        Center.p <- c(Center.p,reads$Center.p)

        names(Left.n)[ids] <- names(reads$Left.n)
        names(Right.p)[ids] <- names(reads$Right.p)
        names(Center.p)[ids] <- names(reads$Center.p)
        names(Center.n)[ids] <- names(reads$Center.n)
      }

      peak.id <- ids[length(ids)]
    }

  }
  if (pairedEnd){
    Reads <-  list(Left.p = Left.p,
                   Left.n = Left.n,
                   Right.p = Right.p,
                   Right.n = Left.n,
                   Center.p = Center.p,
                   Center.n = Center.n,
                   RawTotalCounts = RawTotalCounts,
                   RawCounts.p = RawCounts.p,
                   RawCounts.n = RawCounts.n)
  }else{
    Reads <-  list(Left.p = Left.p,
                   Right.n = Right.n,
                   RawTotalCounts = RawTotalCounts,
                   RawCounts.p = RawCounts.p,
                   RawCounts.n = RawCounts.n)
  }

  return(Reads)
}

sortReads <- function(r, Peaks, pairedEnd = FALSE){
  numPeaks <- length(Peaks)
  Left.p <- vector("list", numPeaks)
  Left.n <- Left.p
  Right.p <- Left.p
  Right.n <- Left.p
  Center.p <- Left.p
  Center.n <- Left.p


  peak.id <- 0
  for (i in 1:length(r)){
    if (length(r[[i]])==0){
      next
    }

    reads <- r[[i]]
    if (length(reads$Left.p)==0&length(reads$Left.n)==0){
      next
    }
    ids <- peak.id+seq(1, length(reads$Left.p))
    Left.p[ids] <- reads$Left.p
    Right.n[ids] <- reads$Right.n

    names(Left.p)[ids] <- names(reads$Left.p)
    names(Right.n)[ids] <- names(reads$Right.n)

    if (pairedEnd){
      Left.n[ids] <- reads$Left.n
      Right.p[ids] <- reads$Right.p
      Center.p[ids] <- reads$Center.p
      Center.n[ids] <- reads$Center.n

      names(Left.n)[ids] <- names(reads$Left.n)
      names(Right.p)[ids] <- names(reads$Right.p)
      names(Center.p)[ids] <- names(reads$Center.p)
      names(Center.n)[ids] <- names(reads$Center.n)
    }
    peak.id <- ids[length(ids)]
  }

  if (pairedEnd){
    Reads <- list(Left.p, Left.n,Right.p, Right.n, Center.p,Center.n)
    names(Reads) <- c('Left.p', 'Left.n','Right.p', 'Right.n',
                      'Center.p','Center.n')

    if (length(Reads$Left.n)!=length(Peaks)){
      stop('some peaks missing')
    }
    if (length(Reads$Right.p)!=length(Peaks)){
      stop('some peaks missing')
    }
    if (length(Reads$Center.p)!=length(Peaks)){
      stop('some peaks missing')
    }
    if (length(Reads$Center.n)!=length(Peaks)){
      stop('some peaks missing')
    }
  }else{
    Reads <- list(Left.p, Right.n)
    names(Reads) <- c('Left.p','Right.n')
  }

  if (length(Reads$Left.p)!=length(Peaks)){
    stop('some peaks missing')
  }

  if (length(Reads$Right.n)!=length(Peaks)){
    stop('some peaks missing')
  }

  ## compare names with peak coordinates
  tmp <- lapply(names(Reads$Left.p),function(n){
    n=unlist(strsplit(n,':'));pos=n[2]})
  st1 <- sapply(tmp,function(pos){
    n=unlist(strsplit(pos,'-'));pos=as.integer(n[1])})
  en1 <- sapply(tmp,function(pos){
    n=unlist(strsplit(pos,'-'));pos=as.integer(n[2])})

  st2 <- start(Peaks)
  en2 <- end(Peaks)
  if (st1 !=st2 || en1!=en2){
    stop("something wrong with peak sorting")
  }


  return(Reads)
}


######################
# getReads: This function uses Rsamtools to collect all short reads on
# chromosome chrom in the bam.file that match to regions defined by
# Peaks. 5' Positions of reads are returned in Reads$p and
# Reads$n for reads mapping to positive or negative strand,
# respectively. Also computes total number of reads mapping to a given peak.
#
#
# INPUT   - Peaks: data frame with cols: chr start stop .  .  strand
#         - chrom: chromosome id
#         - bam.file
#         - pairedEnd
#
# OUTPUT  - Reads$Reads  for single-end reads contains 'Left.p' (Left end of fragments mapping to positive strand)
#                        and 'Right.n' (Right end of fragments mapping to negative strand)
#
#                        for paired-end reads additionally 'Right.p','Left.n', 'Center.p','Center.n',
#
#         - Reads$RawTotalCounts
#         - Reads$RawCounts.p
#         - Reads$RawCounts.n
#
# Gabriele Schweikert
# January 2014
getReads <- function(Peaks, chrom, bam.file, pairedEnd=FALSE){

  header <- scanBamHeader(bam.file)
  if (length(which(names(header[[1]]$targets)==chrom))==0) {
    chrchrom <- paste('chr', chrom, sep='')
    rmchrchrom <- unlist(strsplit(chrom, 'chr'))[2]
    if (length(which(names(header[[1]]$targets)==chrchrom))>0) {
      message('changing ', chrom, ' to ', chrchrom)
      chrom <- chrchrom
    } else if (length(which(names(header[[1]]$targets)==rmchrchrom))>0) {
      message('changing ', chrom, ' to ', rmchrchrom)
      chrom <- rmchrchrom
    } else {
      message('no reads found mapping to ', chrom)
      Names <- vector('list', nrow(Peaks))
      for (i in 1:nrow(Peaks)){
        Names[i] <- paste(Peaks[i,1], ':', Peaks[i,2], '-', Peaks[i,3], sep='')
      }
      Left.p <- vector('list', nrow(Peaks))
      names(Left.p) <- Names
      Left.n <- Left.p
      Right.p <- Left.p
      Right.n <- Left.p
      Center.p <- Left.p
      Center.n <- Left.p
      if  (pairedEnd==TRUE){
        Reads <- list(Left.p, Left.n,Right.p, Right.n,Center.p,Center.n)
        names(Reads) <- c('Left.p', 'Left.n','Right.p', 'Right.n','Center.p','Center.n')
      } else {
        Reads <- list(Left.p, Right.n)
        names(Reads) <- c('Left.p', 'Right.n')
      }
      RawCounts.p <- rep(0,nrow(Peaks))
      RawCounts.n <- rep(0,nrow(Peaks))
      RawTotalCounts <- rep(0,nrow(Peaks))
      Reads <-  list(Reads=Reads,
                     RawTotalCounts=RawTotalCounts,
                     RawCounts.p=RawCounts.p,
                     RawCounts.n=RawCounts.n)
      return(Reads)
    }
  }
  RangesList <- NULL
  IRanges <- NULL
  rm(list=c(RangesList, IRanges))

  which <- RangesList(IRanges(Peaks[,2], Peaks[,3]))
  names(which) <- c(chrom)

  if (pairedEnd==TRUE){

    ## what <- c("rname", "strand", "pos","mpos","isize","seq")
    what <- c("qname","strand", "pos","isize")
    ## isize: This is the TLEN field in SAM Spec v1.4.
    ## Inferred insert size for paired end alignments.

    ## param <- ScanBamParam(which=which, what=what,
    ## scanBamFlag(isProperPair = TRUE))
    ## bam <- scanBam(bam.file, param=param)
    param <- ScanBamParam(which=which, what=what,
                          scanBamFlag(isProperPair = TRUE,
                                      isSecondMateRead = FALSE,
                                      isMinusStrand = FALSE))
    bam <- scanBam(bam.file, param=param)

    if (!all(bam$strand=='+')){
      warning("something wrong with bam filter")
    }
    Left.p <- lapply(bam, function(bam){bam$pos})
    Right.p <- lapply(bam, function(bam){bam$pos + bam$isize})
    Center.p <- lapply(bam, function(bam){bam$pos + .5*bam$isize})

    ## Negative Strand
    param <- ScanBamParam(which=which, what=what,
                          scanBamFlag(isProperPair = TRUE,
                                      isSecondMateRead = TRUE,
                                      isMinusStrand = FALSE))
    bam <- scanBam(bam.file, param=param)
    if (!all(bam$strand=='+')){
      warning("something wrong with bam filter")
    }
    Left.n <- lapply(bam, function(bam){bam$pos})
    Right.n <- lapply(bam, function(bam){bam$pos + bam$isize })
    Center.n <- lapply(bam, function(bam){bam$pos + .5*bam$isize})


    # get total coverage of Peaks
    Counts.p <-  mapply(Left.p, FUN=function(pos){length(pos)})
    Counts.n <-  mapply(Left.n, FUN=function(pos){length(pos)})
    Counts.Right.p <-  mapply(Right.p, FUN=function(pos){length(pos)})
    Counts.Right.n <-  mapply(Right.n, FUN=function(pos){length(pos)})
    if (!identical(Counts.p,Counts.Right.p ) |
        ! identical(Counts.n,Counts.Right.n ))
      stop('')
    RawTotalCounts <- RawCounts.p+RawCounts.n

    Reads <-  list(Left.p = Left.p,
                   Left.n = Left.n,
                   Right.p = Right.p,
                   Right.n = Right.n,
                   Center.p = Center.p,
                   Center.n = Center.n,
                   RawTotalCounts=RawTotalCounts,
                   RawCounts.p=RawCounts.p,
                   RawCounts.n=RawCounts.n)

  } else {
    what <- c("rname", "strand", "pos","seq")
    param <- ScanBamParam(which=which, what=what)
    bam <- scanBam(bam.file, param=param)
    Left.p <- lapply(bam, function(bam){bam$pos[bam$strand=='+']})
    read.Length <- unique(width(bam[[1]]$seq))
    if (length(read.Length)>1){
      warning('reads have different read lengths:')
      print(read.Length)
      read.Length = mean(read.Length)
    }
    Right.n <- lapply(bam, function(bam){bam$pos[bam$strand=='-']+read.Length})

    ## get total coverage of Peaks
    RawCounts.p <-  mapply(Left.p, FUN=function(pos){length(pos)})
    RawCounts.n <-  mapply(Right.n, FUN=function(pos){length(pos)})
    RawTotalCounts <- RawCounts.p+ RawCounts.n

    Reads <-  list(Left.p = Left.p,
                   Right.n = Right.n,
                   RawTotalCounts=RawTotalCounts,
                   RawCounts.p=RawCounts.p,
                   RawCounts.n=RawCounts.n)

  }

  return(Reads)
}
