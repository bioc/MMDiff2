createPeakMatrix <- function(DBA, anchor='start',whichReadPos="Center"){


  P <- DBA$MD$PeakHists[[whichReadPos[1]]]
  nPeaks <-  length(P)
  nSamples <- nrow(P[[1]])
  l <- lengths(P)
  nBins <- max(l)/nSamples

  PeakMatrix <-  vector("list", length(whichReadPos))
  names(PeakMatrix) <- paste('Read', whichReadPos,'-anchoredat-Region',
                             anchor,sep='')

  for (t in 1:length(whichReadPos)){
    P <- DBA$MD$PeakHists[[whichReadPos[t]]]
    H <- vector("list", nSamples)
    names(H) <- DBA$samples$SampleID

    for (i in 1:nSamples){
      Hist <- matrix(NA,nrow=nPeaks,ncol=nBins)
      rownames(Hist) <- names(P)
      colnames(Hist) <- colnames(P[[which.max(l)]])
      for (j in 1:nPeaks){
        if (anchor == 'start'){
          Hist[j,1:ncol(P[[j]])]=P[[j]][i,]
        } else if (anchor == 'end'){
          Hist[j, (nBins-ncol(P[[j]])+1):nBins]=P[[j]][i,]
        } else if (anchor == 'center'){
          Hist[j,(nBins-ncol(P[[j]]))/2+(1:ncol(P[[j]]))]=P[[j]][i,]
        }
      }
      if (anchor == 'end'){
        colnames(Hist) <- as.numeric(colnames(Hist))-
          as.numeric(colnames(Hist)[nBins])
      } else if (anchor == 'center'){
        colnames(Hist) <- as.numeric(colnames(Hist))-
          as.numeric(colnames(Hist)[nBins/2])
      }
      H[[i]] <- Hist
    }
    PeakMatrix[[t]] <- H
  }
  DBA$MD$PeakMatrix <- PeakMatrix
  return(DBA)
}
