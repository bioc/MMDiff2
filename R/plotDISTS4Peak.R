#' plotDISTS4Peak
#'
#' showing all distances for one region
#'
#' @inheritParams DBAmmd-Accessors
#' @inheritParams compDists
#' @inheritParams plotDists
#' @inheritParams plotPeak
#' @param Zoom (DEFAULT: TRUE)
#' @param title (DEFAULT: NULL)
#'
#' @examples
#' dev.off()
#' load(system.file("data/MMD.RData", package="MMDiff2"))
#' plotDISTS4Peak(MMD,Peak.id = '6',dist.method='MMD', whichContrast=1)
#'
#' @export
#
plotDISTS4Peak <- function(MD,Peak.id,dist.method='MMD',whichContrast=1,
                           Zoom=TRUE,xlim=NULL,ylim=NULL,xlog10=TRUE,title=NULL){


  Peak.idx <- match(Peak.id, names(Regions(MD)))

  samples <- Samples(MD)
  Contrast <- MD@Contrasts[[whichContrast]]
  group1 <- names(which(Contrast$group1))
  group2 <- names(which(Contrast$group2))

  DISTs <- Dists(MD,dist.method)
  mCounts <- MD@mCounts

  compNames <- colnames(DISTs)

  within1 <- determineGroupComps(group1,type='within')
  within2 <- determineGroupComps(group2,type='within')
  between <- determineGroupComps(group1,group2,type='between')

  ids1 <- findComps(within1,compNames)
  ids2 <- findComps(within2,compNames)
  ids3 <- findComps(between,compNames)

  W1.y <- as.matrix(DISTs[,ids1])
  W2.y <- as.matrix(DISTs[,ids2])
  B.y <- as.matrix(DISTs[,ids3])

  W1.x <- as.matrix(mCounts[,ids1])
  W2.x <- as.matrix(mCounts[,ids2])
  B.x <- as.matrix(mCounts[,ids3])
  if (xlog10){
    W1.x <- log10(W1.x)
    W2.x <- log10(W2.x)
    B.x <- log10(B.x)
  }

  #### plotting
  if (is.null(xlim)){
    if (Zoom){
      xx <- cbind(B.x,W1.x,W2.x)[Peak.idx,]
      xlim <- c(min(xx[!is.infinite(xx)]),max(xx[!is.infinite(xx)]))
      xlim[1] <- xlim[1] -.2*(xlim[2]-xlim[1])
      xlim[2] <- xlim[2] +.4*(xlim[2]-xlim[1])
    } else {
      xx <- cbind(B.x,W1.x,W2.x)
      xlim <- c(min(xx[!is.infinite(xx)]),max(xx[!is.infinite(xx)]))
    }
  }

  if (is.null(ylim)){
    if (Zoom){
      yy <- cbind(B.y,W1.y,W2.y)[Peak.idx,]
      ylim <- c(min(yy[!is.infinite(yy)]),max(yy[!is.infinite(yy)]))
      ylim[1] <- ylim[1] -.2*(ylim[2]-ylim[1])
      ylim[2] <- ylim[2] +.4*(ylim[2]-ylim[1])
    } else {
      yy <- cbind(B.y,W1.y,W2.y)
      ylim <- c(min(yy[!is.na(yy)]),max(yy[!is.na(yy)]))
    }
  }

  xlab <- paste('mean counts')
  ylab <- paste(dist.method, ' distance between profiles',sep='')

  if (is.null(title)){
    main <- paste('Peak id: ',Peak.id,sep='')
  } else {
    main <- title
  }

  # Add extra space to right of plot area; change clipping to figure
  par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
  if (dist.method=='Pearson'){
    position <- "bottomright"
  }  else if (dist.method=='MMD'){
    position <- "topright"
  } else{
    position <- "topleft"
  }

  plot(B.x[Peak.idx,],B.y[Peak.idx,],xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,
       col = "blue",pch =4, cex = 0.8,cex.main=1.3,cex.lab=1.3,cex.axis=1.3)#,log='x')

  if (is.null(within1)&&is.null(within2)){
    message('no replicates given for within group variance')
  }
  points(W1.x[Peak.idx,],W1.y[Peak.idx,],pch=20,cex=0.8,col='coral')
  points(W2.x[Peak.idx,],W2.y[Peak.idx,],pch=20,cex=0.8,col='orange')


  if (dist.method=='MMD'){
    if (!is.null(Contrast$MMD$fit.Group1)){
      yy=predict(Contrast$MMD$fit.Group1,seq(xlim[1],xlim[2],length.out = 1000))
      points(seq(xlim[1],xlim[2],length.out = 1000),yy,col='coral',t='l',lwd=2)

      d=predict(Contrast$MMD$fit.disp.Group1,seq(xlim[1],xlim[2],length.out = 1000))
      points(seq(xlim[1],xlim[2],length.out = 1000),yy+d,col='coral',t='l',lwd=2,lty=2)
      points(seq(xlim[1],xlim[2],length.out = 1000),yy-d,col='coral',t='l',lwd=2,lty=2)
    }
    if (!is.null(Contrast$MMD$fit.Group2)){
      yy=predict(Contrast$MMD$fit.Group2,seq(xlim[1],xlim[2],length.out = 1000))
      points(seq(xlim[1],xlim[2],length.out = 1000),yy,col='orange',t='l',lwd=2)

      d=predict(Contrast$MMD$fit.disp.Group2,seq(xlim[1],xlim[2],length.out = 1000))
      points(seq(xlim[1],xlim[2],length.out = 1000),yy+d,col='orange',t='l',lwd=2,lty=2)
      points(seq(xlim[1],xlim[2],length.out = 1000),yy-d,col='orange',t='l',lwd=2,lty=2)
    }
  }
  if (Zoom){
    text(B.x[Peak.idx,],B.y[Peak.idx,]-.01,names(ids3),col='blue',cex=.6)
    if (length(ids1)>0){
      text(W1.x[Peak.idx,],W1.y[Peak.idx,]-.01,names(ids1),col='coral',cex=.6)
    }
    if (length(ids2)>0){
      text(W2.x[Peak.idx,],W2.y[Peak.idx,]-.01,names(ids2),col='orange',cex=.6)
    }
  }

  legend(position,inset=c(-.2,0),pch=c(4,20,20),pt.cex=c(1,1,1),
         col=c('blue','coral','orange'),
         c('between groups','within group1', 'within group2'))

  grid()

}
