#' plot Peak
#'
#' This function plots histograms of fragment positions over a pre defined
#' regions of interests / peaks. Can also show occurences of Sequence motifs and
#' annotated objects (e.g. genes).
#'
#' @inheritParams getPeakReads
#' @inheritParams DBAmmd-Accessors
#' @param Peak.id Peak id to specify which Peak to plot.
#' (coresponding to names of Regions(MD))
#' @param Sample.ids which samples to draw. If NULL all samples are drawn.
#' (DEFAULT: NULL)
#' @param NormMethod whether to apply normailzation factors.
#' currently no normalization method implemented (DEFAULT: None)
#' @param plot.input whether to plot input controls (DEFAULT: TRUE)
#' @param Motifs TF binding sites (DEFAULT: NULL)
#' @param Motifcutoff (Default: "80\%")
#' @param anno either a GRanges objects containing annotated objects, e.g. genes,
#' or a list of GRanges Objects. (Default: NULL)
#' @param xaxt (Default: NULL)
#' @param xlim (Default: NULL)
#' @param ylim (Default: NULL)
#'
#' @examples
#' dev.off()
#' data("MMD")
#' plotPeak(MMD,Peak.id='6',plot.input=FALSE)
#'
#' # add annotation (Overlapping genes)
#' data("mm9-Genes")
#' GR <- list(UCSCKnownGenes = GR)
#' plotPeak(MMD, Peak.id='6', plot.input = FALSE, anno=GR)
#'
#' # add TF binding sites
#' library('MotifDb')
#' motifs <- query(query(MotifDb, 'Mmusculus'), 'E2F')
#' plotPeak(MMD, Peak.id='6', plot.input = FALSE,
#'        Motifs=motifs,Motifcutoff="80%")
#'
#' # split peaks by contrast
#' plotPeak(MMD, Peak.id='6', plot.input = FALSE, whichContrast=1,
#'        Motifs=motifs,Motifcutoff="80%",anno=GR)
#'
#'
#' @import BSgenome GenomicRanges Biostrings  S4Vectors
#' @export
#
# TO DO check NormFactors!! from DBA
# check Input is working
# Potentially add seqlogos
# check genome in MD if motif
# make an option to plot strands independently

plotPeak <- function(MD, Peak.id, Sample.ids=NULL,NormMethod=NULL,
                     plot.input=FALSE,whichPos="Center",
                     whichContrast=NULL,
                     Motifs=NULL,
                     Motifcutoff="80%",
                     anno=NULL,
                     xaxt=NULL,xlim=NULL,ylim=NULL){

  if (missing(MD))
    stop("MD is missing")

  Peak.idx <- match(Peak.id, names(Regions(MD)))

  Meta <- metaData(MD)
  PeakBoundary <- Meta$AnaData$PeakBoundary
  hists <- Hists(MD,whichPos)
  samples <- Samples(MD)

  if (is.null(hists)){
    message('whichPos = ', whichPos)
    stop("Peak histograms missing. run 'compHists', first. ")
  }
  if (missing(Peak.id))
    stop("No Peak specified")
  Peak <- hists[Peak.idx]
  PeakName <- names(Peak)
  Peak <- Peak[[1]]

  if (is.null(Sample.ids))
    message("No Samples specified, plotting all samples")

  if (is.null(whichContrast)){
    message("No Contrast specified, plotting all samples in one plot")
    add.contplot=0
  } else add.contplot=1


  if (!is.null(Motifs)){
    if(!class(Motifs)=='MotifList')
      stop("Motifs has to be of class 'MotifList'")
    if (is.null(Genome(MD)))
      stop("genome needed to find motifs")
  }
  if (!is.null(anno)){
    if(!class(anno)=='list' & !class(anno)=='GRanges')
      stop("anno has to be a 'GRanges' object or a list of GRanges objects")
    if (class(anno)=='GRanges')
      anno <- list(anno=anno)
  }

  ################
  ## Define layouts
  ################

  par(oma = c(4,1,0,0) + 0.1,
      mar = c(0,0,1,1) + 0.1)
  leg.size <- .6
  xaxt.1 <- 'n'
  xaxt.2 <- 'n'
  xaxt.3 <- 'n'
  xaxt.4 <- 'n'

  if (is.null(Motifs) & is.null(anno)){

    if (add.contplot==0){
      layout(matrix(c(1,2), 1, 2, byrow = TRUE),
             widths=c(1,3), heights=c(1))
      leg.size <- .4
      xaxt.1 <- xaxt
    } else {
      layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),
             widths=c(1,3), heights=c(1))
      xaxt.2 <- xaxt
    }
  }

  if (xor(is.null(Motifs),is.null(anno))) {
    if (!is.null(Motifs)){
      xaxt.3 <- xaxt
    } else {
      xaxt.4 <- xaxt
    }

    if (add.contplot==0){
      layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE),
             widths=c(1,3), heights=c(2,1))
    } else {
      layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE),
             widths=c(1,3), heights=c(2,2,1))
      leg.size <- 0.8
    }
  }

  if (!is.null(Motifs) & !is.null(anno)){
    xaxt.4 <- xaxt
    if (add.contplot==0){
      layout(matrix(c(1,2,3,4,5,6), 3, 2, byrow = TRUE),
             widths=c(1,3), heights=c(2,1,0.5))
      leg.size <- 0.8
    }else{
      layout(matrix(c(1,2,3,4,5,6,7,8), 4, 2, byrow = TRUE),
             widths=c(1,3), heights=c(2,2,1,0.5))
      leg.size <- 1
    }
  }

  if (is.null(xaxt)){
    xlab=PeakName
  } else{
    xlab=''
  }

  if (is.null(Sample.ids)){
    if (add.contplot==0){
      Sample.ids <- samples$SampleID
      Sample.idx <- 1:length(Sample.ids)
    } else {
      contrast <- MD@Contrasts[[whichContrast]]
      Sample.ids <- names(which(contrast$group1))
      Sample.idx <- which(contrast$group1)

      Sample.ids2 <- names(which(contrast$group2))
      Sample.idx2 <- which(contrast$group2)
      SamplePeaks2 <- Peak[Sample.idx2,,drop=FALSE]
      nSamples2 <- length(Sample.idx2)
    }

  } else{
    Sample.idx <- sapply( Sample.ids,function(id){
      Sample.ix <- which(samples$sampleID==id)})
  }
  SamplePeaks <- Peak[Sample.idx,,drop=FALSE]
  nSamples <- length(Sample.idx)

  ##Normalization
  if (is.null(NormMethod)){
    message('No normalization factors applied')
    Factors <- rep(1,length(Sample.idx))
    if (add.contplot==1){
      Factors2 <- rep(1,length(Sample.idx2))
    }
    ylab <- 'counts'
  } else {
    message('normalization currently not implemented. Set NormMethod=NULL')
    NormFactors <- MD$NormFactors[[NormMethod]][[1]]
    ## check if factor is available for all relevant samples
    if (!all( is.element(rownames(SamplePeaks),names(NormFactors)))){

      tmp <- which(!is.element(rownames(SamplePeaks),names(NormFactors)))
      message('no normalization factor found for sample: ',
              rownames(SamplePeaks)[tmp])
      stop('check normalization factors')
    }
    Factors <- NormFactors[sapply(rownames(SamplePeaks),function(name){
      idx=which(name==names(NormFactors))})]
    ylab='normalized counts'
  }

  ##colors, line type, legend text
  lty <- c(rep(1,nSamples))
  cols <- rainbow(nSamples,start=0,end=0.2)
  leg <- Sample.ids
  title <- ''
  if (add.contplot==1){
    lty2 <- c(rep(1,nSamples2))
    cols2 <- rainbow(nSamples2,start=0.5,end=0.7)
    leg2 <- Sample.ids2
    title <- contrast$name1
    title2 <- contrast$name2
  }


  ##INPUT
  if (plot.input){
    inputs <- unique(samples$ControlID[Sample.idx])
    input.idx <- sapply(inputs,function(id){ix=which(rownames(Peak)==id)})
    if (length(input.idx)==0){
      stop('no inputs given. Check column ControlID
           in your SampleSheet.csv file')
    }
    InputPeaks <- Peak[input.idx,, drop=FALSE]
    leg <- c(leg,inputs)
    cols <- c(cols,rep('black',length(input.idx)))
    lty <- c(lty,rep(2,length(input.idx)))


    if (add.contplot==1){
      inputs2 <- unique(samples$ControlID[Sample.idx2])
      input.idx2 <- sapply(inputs2,function(id){ix=which(rownames(Peak)==id)})
      InputPeaks2 <- Peak[input.idx2,, drop=FALSE]
      leg2 <- c(leg2,inputs2)
      cols2 <- c(cols2,rep('black',length(input.idx)))
      lty2 <- c(lty2,rep(2,length(input.idx)))
    }
    }

  ## get genome coordinates for x-values.

  start <- as.integer(strsplit(strsplit(PeakName,':')[[1]][2], '-')[[1]][1])
  end <- as.integer(strsplit(strsplit(PeakName,':')[[1]][2], '-')[[1]][2])
  chr <- strsplit(PeakName,':')[[1]][1]

  x.coords <- as.integer(colnames(SamplePeaks)) + start - PeakBoundary

  ################
  ### FIND Overlapping MOTIFS
  ################
  if (!is.null(Motifs)) {
    package <- Genome(MD)
    message('checking availability of package ',package)
    available <- suppressMessages(suppressWarnings(
      sapply(package, require, quietly = TRUE, character.only = TRUE,
             warn.conflicts = FALSE)))
    if (available){
      Genome <- getBSgenome(Genome(MD))
      chrom_seq <- Genome[[chr]]
      Seq <- Views(chrom_seq,start- PeakBoundary,end+ PeakBoundary)
      TFBS <- vector("list",length(Motifs))
      names(TFBS) <- names(Motifs)
      TFBS2 <- TFBS
      for (j in 1:length(Motifs)){
        TFBS[[j]] <- matchPWM(Motifs[[j]],Seq,Motifcutoff)
        TFBS2[[j]] <- matchPWM(reverseComplement(Motifs[[j]]),Seq,Motifcutoff)
      }
    } else {
      TFBS = vector("list",0)
      TFBS2 = vector("list",0)
      warning('Genome package is not available, needs to be installed
  first in order to find TF Binding sites')
    }

  }

  ################
  ### FIND Overlapping Annotation
  ################

  if (!is.null(anno)) {
    IRanges <- NULL
    rm(IRanges)

    Reg <- GRanges(seqnames=chr,IRanges(start- PeakBoundary,end+ PeakBoundary))
    GR <- vector("list",length(anno))
    names(GR) <- names(anno)
    for (j in 1:length(anno)){
      ov <- findOverlaps(anno[[j]],Reg)
      GR[[j]] <- anno[[j]][unique(queryHits(ov))]

    }
  }

  if (is.null(ylim)){
    # ylim <- c(min(SamplePeaks)-.4*(max(SamplePeaks)-min(SamplePeaks)),
    # 1.2*max(SamplePeaks))
    if ((add.contplot==0)){
      ylim <- c(min(SamplePeaks),1.2*max(SamplePeaks))
    } else {
      ylim <- c(min(c(SamplePeaks,SamplePeaks2)),
                1.2*max(c(SamplePeaks,SamplePeaks2)))
    }
  }

  if (is.null(xlim)){
    xlim <- c(start- PeakBoundary,end+ PeakBoundary)
  }

  ## first plot legend
  plot.new()

  legend("topleft",legend=leg,col=cols,lty=lty,lwd=1.4,inset=.01,cex=leg.size)

  ### plot Peaks (group1 if contrast given)
  plot(x.coords,SamplePeaks[1,]/Factors[1],
       t='l',col=cols[1],ylim=ylim,xlim=xlim,
       ylab=ylab,main=title,bty='n',
       cex.main=1.5,cex=1,cex.lab=1,cex.axis=1,lwd=1.4,xaxt=xaxt.1)
  if (nSamples>1){
    for (i in 2:nSamples){
      points(x.coords,SamplePeaks[i,]/Factors[i],t='l',col=cols[i],lwd=1.4)
    }
  }
  if (plot.input){
    for (i in 1:length(input.idx)){
      points(x.coords,InputPeaks[i,],lty=2,t='l',col='black',lwd=1.4)
    }
  }

  text(xlim[2]-0.1*(xlim[2]-xlim[1]),
       ylim[2]-0.1*(ylim[2]-ylim[1]),
       paste('Peak id:', Peak.id))
  text(xlim[2]-0.1*(xlim[2]-xlim[1]),
       ylim[2]-0.2*(ylim[2]-ylim[1]),
       paste('Fragment', whichPos))

  grid()

  # plot group2

  if (add.contplot==1){
    ## first plot legend
    plot.new()
    legend("topleft",legend=leg2,col=cols2,lty=lty2,lwd=1.4,inset=.01,
           cex=leg.size)

    ### plot Peaks (group2 if contrast given)
    plot(x.coords,SamplePeaks2[1,]/Factors2[1],
         t='l',col=cols2[1],ylim=ylim,xlim=xlim,
         ylab=ylab,main=title2,bty='n',
         cex.main=1.5,cex=1,cex.lab=1,cex.axis=1,lwd=1.4,xaxt=xaxt.2)
    if (nSamples2>1){
      for (i in 2:nSamples2){
        points(x.coords,SamplePeaks2[i,]/Factors2[i],t='l',col=cols2[i],lwd=1.4)
      }
    }
    if (plot.input){
      for (i in 1:length(input.idx2)){
        points(x.coords,InputPeaks2[i,],lty=2,t='l',col='black',lwd=1.4)
      }
    }
    grid()

  }

  ## ---------
  ## PLOT MOTIFS
  ## ---------
  if (!is.null(Motifs)&& (length(TFBS)>0 | length(TFBS2)>0)){
    #seqLogo(Motifs[[1]])
    ylim=c(0,length(TFBS)+1)
    plot(1, t='l',col='white',ylim=ylim,xlim=c(0,100),
         bty='n',
         cex=1.5,cex.lab=1.5,cex.axis=1.5,yaxt='n',xaxt='n')
    for (j in 1:length(TFBS)){
      text(0,j,names(TFBS[j]),adj=c(0,0),cex=.8)
    }
    plot(1, t='l',col='black',ylim=ylim,xlim=xlim,
         bty='n',
         cex=1.5,cex.lab=1.5,cex.axis=1.5,yaxt='n',xaxt=xaxt.3)

    grid()
    for (j in 1:length(TFBS)){
      tfbs = TFBS[[j]]
      tfbs2 = TFBS2[[j]]
      if (length(tfbs)==0 & length(tfbs2)==0  ) next

      if (!is.null(tfbs) & length(tfbs)>0){
        for (k in 1:length(tfbs)){
          lines(c(start(tfbs)[k],end(tfbs)[k]),c(j,j),lw=5)
        }}

      if (!is.null(tfbs2) & length(tfbs2)>0){
        for (k in 1:length(tfbs2)){
          lines(c(start(tfbs2)[k],end(tfbs2)[k]),c(j,j),lw=5,col='red')
        }
      }
    }
  }

  ####GENE Annotation
  if (!is.null(anno)){
    ylim=c(0,length(GR)+1)
    plot(1, t='l',col='white',ylim=ylim,xlim=c(0,100),
         bty='n',
         cex=1.5,cex.lab=1.5,cex.axis=1.5,yaxt='n',xaxt='n')
    for (j in 1:length(GR)){
      text(0,j,names(GR)[j],adj=c(0,0),cex=1)
    }
    plot(1, t='l',col='black',ylim=ylim,xlim=xlim,
         bty='n',
         cex=1.5,cex.lab=1.5,cex.axis=1.5,yaxt='n',xaxt=xaxt.4)
    grid()
    for (j in 1:length(GR)){
      gr = GR[[j]]

      if (length(gr)==0) next

      for (k in 1:length(gr)){
        lines(c(start(gr)[k],end(gr)[k]),c(j,j),lw=5,col='blue')
        if (as.vector(strand(gr))[k]=='+'){
          points(c(start(gr)[k],end(gr)[k]),c(j,j),cex=3,col='blue',pch='>')
        }
        if (as.vector(strand(gr))[k]=='-'){
          points(c(start(gr)[k],end(gr)[k]),c(j,j),cex=3,col='blue',pch='<')
        }
        if (k==1){
          ov <- intersect(IRanges(start(gr[k]),end(gr[k])),
                          IRanges(xlim[1],xlim[2]))
          text(c(start(ov)+.5*width(ov)),j+0.5,names(gr)[k])
        }

      }
    }
  }
  mtext(text=xlab,side=1,line=3)
}

