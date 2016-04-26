#' plotDists
#'
#' scatterplot showing distances between peaks
#'
#' @inheritParams DBAmmd-Accessors
#' @inheritParams compDists
#' @inheritParams plotPeak
#' @param which.group1 subset samples from group1 (DEFAULT: NULL)
#' @param which.group2 subset samples from group2 (DEFAULT: NULL)
#' @param diff.method which method to use to determine significant peaks
#' (DEFAULT: 'MMD.locfit')
#' @param bUsePval if TRUE p-values instead of FDRs are used (DEFAULT: FALSE)
#' @param th significance threshold for differential called peaks (DEFAULT: 0.1)
#' @param title an overall title for the plot (DEFAULT: NULL)
#' @param what which dists to overlay: 1: only between group distances,
#' 2: between and within group distances, 3: between and within group distances,
#' and significant peaks highlightend   (DEFAULT: 3)
#' @param xlim specify x range (DEFAULT: NULL)
#' @param ylim specify y range (DEFAULT: NULL)
#' @param xlog10 should x range be plotted in log10 scale (DEFAULT: TRUE)
#' @param Peak.IDs Highlight specific subset of peaks (DEFAULT: NULL)
#' @param withLegend (DEFAULT: TRUE)
#' @param shiny_df_opt Option returns a dataframe for shiny (DEFAULT: FALSE)
#'
#' @examples
#' data("MMD")
#' plotDists(MMD, whichContrast=1)
#'
#' @export
#
plotDists <- function(MD,
                      dist.method='MMD',
                      whichContrast=1,
                      which.group1=NULL,
                      which.group2=NULL,
                      diff.method='MMD.locfit',
                      bUsePval=FALSE,
                      th=0.1,
                      title=NULL,
                      what=3,
                      xlim=NULL,
                      ylim=NULL,
                      xlog10=TRUE,
                      Peak.IDs=NULL,
                      withLegend=TRUE,
                      shiny_df_opt=FALSE){

  samples <- Samples(MD)
  Contrast <- MD@Contrasts[[whichContrast]]

  DISTs <- Dists(MD,dist.method)
  mCounts <- MD@mCounts


  if (is.element(diff.method,names(Contrast))){
    M <- Contrast[[diff.method]]
    if (diff.method =='MMD.locfit'){
      if (bUsePval){
        Pvals <- M$mmd.Table[,'pval']
      } else {
        Pvals <- M$mmd.Table[,'padj']
      }
      Idx <- which(Pvals<th)
      i.o <- which(M$outlier==TRUE)
      Idx <- setdiff(Idx,i.o)
    } else if (diff.method =='edgeR'){
      # rep <- dba.report(DBA,contrast,method=diff.method,bUsePval=bUsePval)
      # Idx <- as.numeric(names(rep) )
    }
  } else Pvals <- NULL

  compNames <- colnames(DISTs)

  group1 <- names(which(Contrast$group1))
  group2 <- names(which(Contrast$group2))
  if (!is.null(which.group1)){
    group1 <- intersect(group1,which.group1)
  }
  if (!is.null(which.group2)){
    group2 <- intersect(group2,which.group2)
  }



  within1 <- determineGroupComps(group1,type='within')
  within2 <- determineGroupComps(group2,type='within')
  between <- determineGroupComps(group1,group2,type='between')

  ids1 <- findComps(within1,compNames)
  ids2 <- findComps(within2,compNames)
  ids3 <- findComps(between,compNames)

  W1.y <- as.matrix(DISTs[,ids1])
  W2.y <- as.matrix(DISTs[,ids2])
  W.y <- rowMeans(cbind(W1.y,W2.y))
  B.y <- rowMeans(as.matrix(DISTs[,ids3]))
  #B.y <- DISTs[,ids3[1]]
  #W.y <- W1.y[,1]

  W1.x <- as.matrix(mCounts[,ids1])
  W2.x <- as.matrix(mCounts[,ids2])

  W.x <- rowMeans(cbind(W1.x,W2.x))
  B.x <- rowMeans(as.matrix(mCounts[,ids3]))
  if (xlog10){
    W.x <- log10(W.x)
    B.x <- log10(B.x)
  }
  #B.x <- mCounts[,ids3[1]]
  #W.x <- W1.x[,1]


  #### plotting
  if (is.null(xlim)){
    xx <- cbind(B.x,W.x)
    xlim <- c(min(xx[!is.infinite(xx)]),max(xx[!is.infinite(xx)]))
  }
  if (is.null(ylim)){
    yy <- cbind(B.y,W.y)
    ylim <- c(min(yy[!is.na(yy)]),max(yy[!is.na(yy)]))
  }


  xlab <- paste('log10 mean counts')
  ylab <- paste(dist.method, ' distance between profiles',sep='')

  if (is.null(title)){
    main <- paste(Contrast$name1,' vs ',Contrast$name2,
                  ' (nPeaks = ',dim(mCounts)[1],')',sep='')
  } else {
    main <- title
  }

  if (!is.null(Pvals)){
    if (bUsePval){
      sigred='pval'
    } else sigred='FDR'
    leg3 <- paste(sigred,'<', th, ': n = ',length(Idx),
                  ' (',diff.method,')',
                  sep='')

    # leg3 <- paste(sigred,'<', thresh, ' (',colnames(Pvals)[idx],'):
    # n = ',length(Idx),sep='')
  }else{
    leg3 <- ''
  }

  if (dist.method=='Pearson'){
    position <- "bottomright"
  }  else if (dist.method=='MMD'){
    position <- "topright"
  } else {
    position <- "topright"
  }

  plot(B.x,B.y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
       col = "blue",pch = 4, cex = 0.6,cex.main=1.3,cex.lab=1.3,cex.axis=1.3)
  # smoothScatter(Ns,rowMeans(B),
  # xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,col = "blue")
  between_df <- data.frame(cbind(B.x, B.y,
                                 rep(0, length.out=length(B.x)),
                                 rep(0, length.out=length(B.x))),
                           stringsAsFactors = FALSE)
  colnames(between_df) <- c("means","distance","condition","sig")

  if (is.null(within1)&&is.null(within2)){
    message('no replicates given for within group variance')
  } else if (what>1){
    points(W.x,W.y,pch=20,cex=0.6)
    within_df <- data.frame(cbind(W.x, W.y,
                                  rep(1,length.out=length(W.x)),
                                  rep(0, length.out=length(W.x))),
                            stringsAsFactors = FALSE)
    colnames(within_df) <- c("means","distance","condition","sig")
  }
  if (!is.null(Pvals)&&length(Idx)>0 & what>2){
    points(B.x[Idx],B.y[Idx],col='red',pch=21,cex=0.6)
    between_df$sig[Idx] <- 1
  }



  if (diff.method=='MMD.locfit'){
    N <- seq(max(0,min(c(B.x,W.x))),max(c(B.x,W.x)),length.out = 1000)
    log10N <- N
    if (!xlog10){
      log10N <- log10(N)
    }
    if (!is.null(M$fit.Group1)){

      #data <- data.frame(N=10^(seq(xlim[1],xlim[2],length.out = 1000)))
      yy <- predict(M$fit.Group1,log10N)
      points(N,yy,col='orange',t='l',lwd=2)

      yyd <- predict(M$fit.disp.Group1,log10N)
      points(N,yy+yyd,col='yellow',t='l',lwd=2,lty=2)
      points(N,yy-yyd,col='yellow',t='l',lwd=2,lty=2)
    }
    if (!is.null(M$fit.Group2)){
      yy=predict(M$fit.Group2,log10N)
      points(N,yy,col='orange',t='l',lwd=2)

      yyd <- predict(M$fit.disp.Group2,log10N)
      points(N,yy+yyd,col='yellow',t='l',lwd=2,lty=2)
      points(N,yy-yyd,col='yellow',t='l',lwd=2,lty=2)

    }
  }


  if (!is.null(Peak.IDs)){
    Peak.idxs <- match(Peak.IDs, names(Regions(MD)))
    points(B.x[Peak.idxs],B.y[Peak.idxs],col='green',pch=21,cex=1)
  }
  if (withLegend){
    legend(position,pch=c(4,20,21),pt.cex=c(1,1,1),
           col=c('blue','black','red'),
           c('between groups','within groups',leg3),
           inset=.01)
  }
  grid()

  if (shiny_df_opt == TRUE) {
    rownames(within_df) <- paste('W-',names(Regions(MD)),sep='')
    rownames(between_df) <- paste('B-',names(Regions(MD)),sep='')
    out_df <- rbind(within_df, between_df)
    colnames(out_df) <- c("means","distance","condition","sig")
    return(out_df)
  }
  #legend(position,fill=c('blue','white','white'),
  #border=c('black','white','white'),pch=c(20,4,21),pt.cex=c(0.2,0.6,0.8),
  #col=c('black','black','red'),c('between groups','within groups',
  #paste('p<0.05 (',colnames(Pvals)[idx],' )')))


}
