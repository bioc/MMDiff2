# compHist
# addHists
# compPeakMatrix
# determineGroupComps
# findComps

compHist <- function(Pos, PeakLengths, bin.length=20){

  ## N <- names(Pos)
  ## aa <- strsplit(N,':')
  ## S <- sapply(aa, function(a){strsplit(a[2], '-')})
  #breaks <- lapply(PeakLengths, function(L){seq(from=1, to=L, by=bin.length)})
  breaks <- vector('list',length=length(PeakLengths))
  for (i in 1:length(PeakLengths)){
    breaks[[i]]<- seq(from=1, to=PeakLengths[i], by=bin.length)
  }
  Counts <- vector("list", length(Pos))
  names(Counts) <- names(Pos)
  Counts <- Counts
  Mids <- Counts

  for (i in 1:length(Pos)){ #loop over peaks
    if (length(Pos[[i]])==0){
      Counts[[i]] <- rep(0, length(breaks[[i]])-1)
      Mids[[i]] <- hist(breaks[[i]], breaks[[i]], plot=FALSE)$mids
    } else {
      P <- Pos[[i]][Pos[[i]]>min(breaks[[i]])&Pos[[i]]<max(breaks[[i]])]
      h <- hist(P, breaks=breaks[[i]], plot = FALSE)
      Counts[[i]] <- h$counts
      Mids[[i]] <- h$mids
    }
  }

  Hists <- list(Counts, Mids)
  names(Hists) <- c('Counts', 'Mids')
  return(Hists)
}

addHists <- function(Hist1,Hist2){
  if (!identical(names(Hist1),  names(Hist2))){
    stop("trying to add different histograms")}
  Hist <- vector("list", length(Hist1))
  names(Hist) <- names(Hist1)
  for (i in 1:length(Hist1)){
    Hist[[i]] <- Hist1[[i]]+Hist2[[i]]
  }
  return(Hist)
}



determineGroupComps <- function(group1,group2,type){

  n1 <- length(group1)

  if (identical(type,'within')){
    if (n1<2){
      within <- NULL
      return(within)
    }

    ID <- which(upper.tri(matrix(1,n1,n1), diag = FALSE),arr.ind=TRUE)
    CompIDs <- matrix(0,ncol=nrow(ID),nrow=ncol(ID))
    CompIDs[1,] <- group1[ID[,1]]
    CompIDs[2,] <- group1[ID[,2]]

    within <- paste(CompIDs[1,],CompIDs[2,],sep=' vs ')
    return(within)
  } else {
    n2 <- length(group2)


    st <- 0
    nComps <- n1*n2
    CompIDs <- matrix(0,2,nComps)
    for (i in 1:n1){
      CompIDs[1,st+(1:n2)] <- group1[i]
      CompIDs[2,st+(1:n2)] <- group2

      st <- st+n2

    }
    between <- paste(CompIDs[1,],CompIDs[2,],sep=' vs ')
    return(between)

  }
}

findComps <- function(Names,compNames,method='MMD'){

  ids <- sapply(Names,function(id){ids=which(compNames==id)})
  i <- which(mapply(length,ids)==0)
  #try to swap comparison
  if (length(i)>0){
    N <- Names[i]
    temp <- strsplit(N,' vs ')
    Nrev <- sapply(temp,function(t){paste(t[2],t[1],sep=' vs ')})
    idsrev <- sapply(Nrev,function(id){ids=which(compNames==id)})
    ids[i] <-  idsrev
  }
  i <- which(mapply(length,ids)==0)
  if (length(i)>0){
    warning(method, ' distances for ',Names[i],' not found')
  }
  ids <- unlist(ids)
  return(ids)
}

######################
# compPeakMatrix takes the Histograms in fieldName='Counts' of Hists
# and reorganises the data such that for each peak a matrix is
# generated containing histograms from all samples
#
# INPUT   - Hists[[fieldName]]
#         - fieldName='Counts'
#
# OUTPUT  - PeakHists
#
#
# Gabriele Schweikert
# August 2012

compPeakMatrix <- function(Hists, fieldName='Counts'){

  nPeaks <- lapply(Hists, function(h){length(h[[fieldName]])})
  nPeaks <- as.integer(unique(nPeaks))
  if (length(nPeaks)!=1){
    stop('different peak sets per sample!')
  }


  nSamples <- length(Hists)
  Samples <- names(Hists)

  PeakHists <- vector("list", length(Hists[[1]][[fieldName]]))
  names(PeakHists) <-  names(Hists[[1]][[fieldName]])
  for (i in 1:nPeaks){
    Mids <- Hists[[1]]$Mids[i]
    Counts <- matrix(0, nrow=nSamples, ncol=length(Mids[[1]]))
    colnames(Counts) <- Mids[[1]]
    rownames(Counts) <- Samples
    for (j in 1:nSamples){
      Counts[j,] <- Hists[[j]][[fieldName]][[i]]
      if (!identical(Mids[[1]], Hists[[j]]$Mids[i][[1]])){
        stop('Mid points do not agree in different samples')
      }

    }
    PeakHists[[i]] <- Counts
  }

  return(PeakHists)

}


