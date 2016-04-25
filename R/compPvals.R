#' compute p-values
#'
#' This function determines peak-specific p-values based on distances between
#' sample histograms.
#'
#' @inheritParams compDists
#' @param diff.method (DEFAULT: 'MMD.locfit')
#' @return DBAmmd object with updated Contrast Dists
#'
#' @examples
#'
#' ## Example using a small data set provided with this package:
#'
#' data("MMD")
#' MMD.1 <- compPvals(MMD)
#'
#' @seealso \code{\link{DBAmmd}}, \code{\link{plotDists}},\code{\link{compDists}}
#'
#'
#' @import Biobase locfit
#' @importFrom stats p.adjust residuals predict pchisq pt
#' @export
#
#
# subfunctions:
#         - compute.pvalues
#         - determineGroupComps
#         - findComps
#
# Gabriele Schweikert
# March 2016

compPvals <- function(MD, dist.method='MMD',diff.method='MMD.locfit'){

  if (missing(MD))
    stop("DBAmmd is missing")


  Meta <- metaData(MD)

  DISTs <- Dists(MD,dist.method)
  mCounts <- MD@mCounts
  contrasts <- MD@Contrasts
  if (is.null(contrasts)) {
    stop('contrasts missing')
  }


  for (c in 1:length(contrasts)){
    contrast <- contrasts[[c]]
    group1 <- names(which(contrasts[[c]]$group1))
    group2 <- names(which(contrasts[[c]]$group2))

    compNames <- colnames(DISTs)

    within1 <- determineGroupComps(group1,type='within')
    within2 <- determineGroupComps(group2,type='within')
    between <- determineGroupComps(group1,group2,type='between')

    if (is.null(within1)&&is.null(within2)){
      stop('no replicates given to estimate within group variance')}
    if (is.null(within1)){
      warning('no replicates given for group1')}
    if (is.null(within2)){
      warning('no replicates given for group2')}

    group1.ids <- findComps(within1,compNames,dist.method)
    group2.ids <- findComps(within2,compNames,dist.method)
    bw.group.ids <- findComps(between,compNames,dist.method)

    message("Computing p-values \n")

    if (diff.method=='MMD.locfit'){
      P <- compute.pvalues(mCounts, MMD=DISTs, group1.ids, group2.ids,
                           bw.group.ids)
      contrast$MMD.locfit <- list(fit.Group1=P$fit.Group1, fit.Group2=P$fit.Group2,
                           fit.disp.Group1=P$fit.disp.Group1,
                           fit.disp.Group2=P$fit.disp.Group2,
                           group1=P$group1,
                           group2=P$group2,
                           between.gr1=P$between.gr1,
                           between.gr2=P$between.gr2,
                           mmd.Table=P$mmd.Table,
                           outlier=P$outlier)
    }
    contrasts[[c]] <- contrast
  }
  MD@Contrasts=contrasts
  return(MD)
}



compute.pvalues <- function(Counts, MMD, group1.ids, group2.ids, bw.group.ids ){


  Counts <- log10(Counts+1)
  quants <- quantile(as.vector(Counts),seq(0,1,0.01))

  mmd.Table = matrix(1,dim(Counts),3)
  colnames(mmd.Table) = c('id','pval','padj')
  mmd.Table[,'id'] =seq(1,dim(Counts)[1])
  # FILTER
  idx = which(rowMin(Counts) > quants[2] &
                rowMax(Counts) < quants[length(quants)-1] &
                ! is.na(rowSums(MMD)))
  C = Counts[idx,]
  M = MMD[idx,]


  if (!is.null(group1.ids)){
    # Regression MMD on Counts
    data = data.frame(N = rowMeans(as.matrix(C[,group1.ids])),
                      D = rowMeans(as.matrix(M[, group1.ids])))


    fit.Group1 <- locfit.robust(D ~ N, data=data)
    disp <- abs(data$D - predict(fit.Group1,data$N))
    #fit.Group1 <- lm(D ~ 1/N, data=data)
    #disp <- abs(data$D - predict.lm(fit.Group1,data))
    data <- cbind(data,disp = disp)
    fit.disp.Group1 <- locfit.robust(disp ~ N, data=data)
    #fit.disp.Group1 <- lm(disp ~ 1/(10^N), data=data)
    #yy=predict.lm(fit.disp.Group1,data)
    ## find outlier
    group1 = detpvals(MMD,Counts,fit.Group1,fit.disp.Group1,group1.ids)
    Group1 = group1$combined

    between.gr1 <- detpvals(MMD,Counts,fit.Group1,fit.disp.Group1,bw.group.ids)

  }  else {
    fit.Group1 <- NULL
    fit.disp.Group1 <- NULL
    between.gr1 <- NULL
    Group1 = matrix(nrow=nrow(Counts),ncol=0)
  }

  if (!is.null(group2.ids)){
    data = data.frame(N = rowMeans(as.matrix(C[,group2.ids])),
                      D = rowMeans(as.matrix(M[, group2.ids])))

    fit.Group2 <- locfit.robust(D ~ N, data=data)
    disp <- abs(data$D - predict(fit.Group2,  data$N))
    #fit.Group2 <- lm(D ~ 1/N, data=data)
    #disp <- abs(data$D - predict.lm(fit.Group2,  data))
    data <- cbind(data,disp = disp)
    fit.disp.Group2 <- locfit.robust(disp ~ N, data=data)

    group2 = detpvals(MMD,Counts,fit.Group2,fit.disp.Group2,group2.ids)
    Group2 = group2$combined
    between.gr2 <- detpvals(MMD,Counts,fit.Group2,fit.disp.Group2,bw.group.ids)
    if (!is.null(group1.ids)){
      comb <- cbind(between.gr1$combined,between.gr2$combined)
      comb[is.na(comb)] <- 10
      comb <- rowMax(comb)
    } else {
      comb <- between.gr2$combined
    }

  } else {
    fit.Group2 <- NULL
    fit.disp.Group2 <- NULL
    group2 = NULL
    Group2 <- matrix(nrow=nrow(Counts),ncol=0)
    between.gr2 <- NULL
    comb <- between.gr1$combined
    comb[is.na(comb)] <- 10
  }
  ##   p-vals


  fdr <- comb
  fdr[fdr!=10] <- p.adjust(fdr[fdr!=10] ,method="BH")
  mmd.Table[,'pval'] <- comb
  mmd.Table[,'padj'] <- fdr

  withinGroup.pvals = cbind(p.adjust(Group1,method="BH"),
                            p.adjust(Group2,method="BH"))
  withinGroup.pvals[is.na(withinGroup.pvals)]=10
  w.p  = rowMin(withinGroup.pvals)
  pvals = mmd.Table[,'pval']
  outlier =  rep(FALSE,length(w.p))
  outlier[pvals>w.p & w.p<0.1]=TRUE


  #I=sort.int(mmd.Table[,'padj'] ,
  # na.last = NA, decreasing = FALSE, index.return = TRUE)
  # mmd.Table = mmd.Table[I$ix,]

  P=list(fit.Group1=fit.Group1, fit.Group2=fit.Group2,
         fit.disp.Group1=fit.disp.Group1, fit.disp.Group2=fit.disp.Group2,
         mmd.Table=mmd.Table,
         group1=group1,
         group2=group2,
         between.gr1=between.gr1,
         between.gr2=between.gr2,
         outlier=outlier)
  return(P)
}

detpvals <- function(M,C,fit,fit.disp,which.comps){
  p <- matrix(NA,dim(M)[1],length(which.comps))
  s <- summary(fit)
  #N <- length(s$residuals)
  #K <- s$df[2]
  N <- summary(fit)$n
  K <- fit$dp['df2']
  RSDR <- quantile(abs(residuals(fit)),.6827)*N/(N-K)
  for (i in 1:length(which.comps)) {
    pred.y <- predict(fit,C[,which.comps[i]])
    #pred.y <- predict.lm(fit,data.frame(N=C[,which.comps[i]]))
    Dist <- M[,which.comps[i]]-pred.y
    DR <- predict(fit.disp,C[,which.comps[i]])
    t1 <- Dist/RSDR
    t2 <- Dist/DR
    p[,i] <- pt(t2,df=N-K,lower.tail = FALSE)
  }
  W <- -2*rowSums(log(p))
  pp <- 1-pchisq(W,2*length(which.comps))
  #Dist = rowMeans(Dist)
  #t <- Dist/RSDR
  #p <- pt(t,df=N-K,lower.tail = FALSE)
  P <- list(uncombined=p,combined=pp)
  return(P)
}
