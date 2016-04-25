#' Compute distances between Peaks
#'
#' This function computes pairwise distances between histograms
#' according to the dist.method (MMD, KS). For large data sets it is a bit
#' time consuming.
#'
#' @inheritParams DBAmmd-Accessors
#' @inheritParams getPeakReads
#' @param sigma  sigma parameter of the RBF kernel,
#' determining the distance (along the genome) at which fragment counts
#' decorrelate. If set to NULL, 100 random Peaks are used to determine sigma
#' heuristically according to the method described in the MMDiff paper [1].
#' (DEFUALT: NULL)
#'
#' @return DBAmmd object with updated slot Dists
#' @seealso \code{\link{DBAmmd}}, \code{\link{plotDists}},
#' \code{\link{plotDISTS4Peak}}, \code{\link{compPvals}}
#'
#' @examples
#'
#' ## Example using a small data set provided with this package:
#'
#' data("MMD")
#' MMD.1 <- compDists(MMD)
#'
#' # To inspect the computed distances:
#' D <- Dists(MMD.1,dist.method='MMD')
#' head(D)
#'
#' # To analyse the result:
#' plotDists(MMD.1)
#'
#' @importFrom stats quantile median ks.test
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#' @author Gabriele Schweikert \email{G.Schweikert@ed.ac.uk}
#' @references [1] Schweikert et al.  BMC Genomics 2013
#' ...
# subfunctions:
#         - compDist.R
#         - mmd.R
#
# check parameters
#
# Gabriele Schweikert
# March 2016




compDists <- function(MD, dist.method='MMD',sigma=NULL,
                      run.parallel=TRUE){

  message('checking parameters...')

  if (missing(MD))
    stop("DBAmmd object")

  Meta <- metaData(MD)
  PeakBoundary <- Meta$AnaData$PeakBoundary
  Peaks <- Regions(MD)

  if (length(Peaks)<200){
    numPeaksForSigma <- length(Peaks)/2-1
  } else
    numPeaksForSigma <- 100

  #R <- Reads(MD)
  #if (!is.element('Reads',names(MD))){
  #    stop('No reads provided, call getPeakReads first')
  #}
  #if (!is.element('Pos.Center',names(MD$Reads))){
  #    stop('Read Centers not provided, call something first')
  #}

  Counts <- Counts(MD)
  Pos.C <- Reads(MD,'Center')

  # flanking regions
  Flanks <- list(ncol(Counts))
  for (i in 1:ncol(Counts)){
    N <- median(Counts[,i]/width(Peaks)*PeakBoundary)
    Flanks[[i]] <- sample(PeakBoundary,N,replace=TRUE)
  }

  ## ----------------
  ## establish all pairwise comparisons
  ## -----------------

  nSamples <- numSamples(MD)
  sampleIDs <- Samples(MD)$SampleID

  nComps <- (nSamples^2-nSamples)/2
  CompIDs <- matrix(0,2,nComps)
  CompIDxs <- matrix(0,2,nComps)
  st <- 0
  for (j in 1:(nSamples-1)){
    CompIDs[1,st+1:(nSamples-j)] <- sampleIDs[j]
    CompIDs[2,st+1:(nSamples-j)] <- sampleIDs[(j+1):nSamples]

    CompIDxs[1,st+1:(nSamples-j)] <- j
    CompIDxs[2,st+1:(nSamples-j)] <- (j+1):nSamples

    st <- st+nSamples-j
  }

  if (dist.method=='MMD'){
    ## ----------------
    ## 1. estimate sigma
    ## ----------------
    message('estimating sigma...')

    Sigma <- matrix(0,numPeaksForSigma,nSamples )
    colnames(Sigma) <- sampleIDs

    for (i in 1:nSamples){
      C <- Counts[,i]
      peak.ids <- which(C>quantile(C)[2]&C<quantile(C)[4])
      peak.ids <- sample(peak.ids, numPeaksForSigma)
      R <- Pos.C[[i]][peak.ids]
      for(j in 1:numPeaksForSigma){
        Sigma[j,i] <- chooseSigma(R[[j]],R[[j]])
      }
    }
    summary(Sigma)
    sigma <- median(as.vector(Sigma))

    ## ----------------
    ## 2. precompute Kernel matrix
    ## ----------------

    message('pre-computing Kernel matrix...')
    L <- max(width(Peaks))+2*PeakBoundary
    KernelMatrix <- compKernelMatrix(seq(1,L),sigma)

  } #end if (dist.method=='MMD')

  ## ----------------
  ## 3. compute DISTs for all peaks and all sample pairs
  ## ----------------


  message(paste('computing', nComps ,'pair-wise distances...'))
  CompNames <- paste(CompIDs[1,],CompIDs[2,],sep=' vs ')
  D <- matrix(NA,length(Peaks),nComps)
  MCounts <- matrix(NA,length(Peaks),nComps)

  set.seed(1236436)
  for (i in 1:nComps){
    i1 <- CompIDxs[1,i]
    i2 <- CompIDxs[2,i]
    message(paste('computing distances for',
                  CompIDs[1,i],'vs',  CompIDs[2,i]))
    #MCounts[,i] <- rowMeans(Counts[,c(i1,i2)])
    MCounts[,i] <- rowMax(Counts[,c(i1,i2)])
    pb <- txtProgressBar(min=1,max=length(Peaks),style=3)
    Ls <- width(Peaks)+PeakBoundary

    for (j in 1:length(Peaks)){
      posA <- Pos.C[[i1]][[j]]
      posB <- Pos.C[[i2]][[j]]

      if (dist.method=='MMD'){
        if (length(posA)>2 | length(posB)>2){

          posA <- c(posA,Flanks[[i1]],Ls[j]+Flanks[[i1]])
          posB <- c(posB,Flanks[[i2]],Ls[j]+Flanks[[i2]])

          D[j,i] <- mmd(posA,posB,KernelMatrix)
        }
      } else if (dist.method=='KS'){
        KS <- ks.test(posA,posB)
        D[j,i] <- KS$statistic[[1]]
      }


      rm(posA,posB)

      setTxtProgressBar(pb,j)
    }  # end loop over Peaks

  } # end loop over comps

  #D.L[is.na(D.L)] = D.R[is.na(D.L)]
  #D.R[is.na(D.R)] = D.L[is.na(D.R)]
  #DISTs <- (D.L+D.R)/2
  DISTs <- D
  CompNames <- paste(CompIDs[1,],CompIDs[2,],sep=' vs ')
  colnames(DISTs) <- CompNames
  rownames(DISTs) <- names(Pos.C[[1]])
  colnames(MCounts) <- CompNames
  rownames(MCounts) <- names(Pos.C[[1]])
  MD@mCounts <- MCounts
  if (dist.method=='MMD'){
    Meta$AnaData$Flanks <- Flanks
    Meta$AnaData$MMDKernelSigma <- Sigma
    Meta$AnaData$MMDKernelFinalSigma <- sigma
    Meta$AnaData$KernelMatrix <- KernelMatrix
  }
  ## ----------------
  ## 4. preparing Output
  ## ----------------

  message('preparing Output...')


  MD@MetaData <- Meta
  MD <- setDists(MD,DISTs,dist.method)

  return(MD)
}



chooseSigma <- function(x,y){
  x=as.matrix(x,nrow=(length(x)))
  y=as.matrix(y,nrow=(length(y)))
  Kxy=vectorized_pdist.sq(x,y)

  mdist <- median(Kxy[Kxy!=0])

  mdist <- sqrt(mdist/2)
  sigma <- mdist

  # X <- seq(1,max(c(x,y)))
  # D <- as.matrix(dist(X,diag=TRUE,upper=TRUE))
  # DISTs <- D[x,y]
  # Kxy <- DISTs^2

  if (sigma ==0){
    sigma <- 1}
  return(sigma)
}


compKernelMatrix <- function(x,sigma){
  ## should do the same as rbf_dot
  #DISTs <- as.matrix(dist(x),diag = TRUE, upper = TRUE)
  #K <-  exp(-1/2/sigma^2 * DISTs^2)
  x=as.matrix(x,nrow=(length(x)))
  DISTs.sq <- vectorized_pdist.sq(x,x)
  K <-  exp(-1/2/sigma^2 * DISTs.sq)
  return(K)
}



mmd <- function(posA,posB,KernelMatrix) {

  m <- length(posA)
  n <- length(posB)

  Kxx <- KernelMatrix[posA,posA]
  Kyy <- KernelMatrix[posB,posB]
  Kxy <- KernelMatrix[posA,posB]


  ## Get test statistic (unbiased)
  # MMDf = ( Kxx + Kyy - Kxy - t(Kxy) )
  # MMDf = MMDf - diag(diag(MMDf))
  # testStat = 1/m/(m-1) * sum(sum( MMDf ));
  # testStat = testStat * m; #null distirbution on m*MMD

  ## MMD statistic. Here we use biased
  ## v-statistic. NOTE: this is m * MMD_b
  # testStat = 1/m * sum(sum(Kxx + Kyy - Kxy - t(Kxy)))
  biased <- sqrt(sum(Kxx)/(m*m) +  sum(Kyy)/(n * n) - 2/m/n * sum(Kxy))
  # diag(Kxx) <- 0
  # diag(Kxy) <- 0
  # unbiased <- sqrt(sum(Kxx)/(m*(m-1)) +  sum(Kyy)/(n * (n-1)) - 2/m/n * sum(Kxy))
  return(biased)
}

vectorized_pdist.sq <- function(A,B) {
  #from Alex Smola
  an = apply(A, 1, function(rvec) crossprod(rvec,rvec))
  bn = apply(B, 1, function(rvec) crossprod(rvec,rvec))

  m = nrow(A)
  n = nrow(B)

  tmp = matrix(rep(an, n), nrow=m)
  tmp = tmp +  matrix(rep(bn, m), nrow=m, byrow=TRUE)
  d.sq = tmp - 2 * tcrossprod(A,B)
  return(d.sq)
}



