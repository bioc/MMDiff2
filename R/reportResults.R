#' report results
#'
#' retrieve results of differential binding analysis
#'
#'
#' @inheritParams compPvals
#' @inheritParams plotDists
#' @param rm.oulier if TRUE, significant peaks with high
#' within-group distances are not reported. (DEFAULT: TRUE)
#' @examples
#'
#' data("MMD")
#' res <- reportResults(MMD)
#'
#' @export

reportResults <- function(MD, diff.method='MMD.locfit', th=0.1,
                          whichContrast=1,rm.oulier=TRUE,bUsePval=FALSE){

  c <- Contrast(MD)[[whichContrast]]
  c.m <- c[[diff.method]]

  if (diff.method=='MMD.locfit'){
    if (!bUsePval){
      idx.sig <- which(c.m$mmd.Table[,'padj'] < th)
    } else{
      idx.sig <- which(c.m$mmd.Table[,'pval'] < th)
    }
    if (rm.oulier){
      i.o <- which(c.m$outlier==TRUE)
      idx.sig <- setdiff(idx.sig,i.o)
    }
    res <- Regions(MD)[idx.sig]
    tab <- c.m$mmd.Table[idx.sig,]
    ii <- sort(tab[,'padj'],index.return = TRUE)
    res <- res[ii$ix]
    tab <- tab[ii$ix,]
    res$padj <- tab[,'padj']
    res$pval <- tab[,'padj']
  }
  return(res)
}