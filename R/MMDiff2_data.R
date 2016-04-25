#' Peaks for Cfp1-data set
#'
#' Subset of MACS called Peaks for Cfp-1 data set. Consensus Peaks were created
#' using diffBind (see below).
#'
#' @name Cfp1-Peaks
#' @docType data
#' @references data taken from Clouaire et al., Genes and Development, 2012.
#' @usage data('Cfp1-Peaks')
#' @keywords data
#' @format contains Peaks, a GRanges object with 500 ranges and
#' 3 metadata columns
#' @examples
#' # data was created as follows:
#' \dontrun{
#' library('MMDiffBamSubset')
#' dataDir <- system.file("extdata", package="MMDiffBamSubset")
###  DiffBind
#' library('DiffBind')
#' olddir <- setwd(dataDir)
#' DBA <- dba(sampleSheet="Cfp1.csv", minOverlap=3)
#' Peaks <- dba.peakset(DBA, bRetrieve = TRUE)
#' DBA <- dba.count(DBA, minOverlap=3)
#' setwd(olddir)
#' peaks <- dba.peakset(DBA, bRetrieve=TRUE)
#' C <- Counts(MMD)
#' idx <- which(C[,1]>150 & C[,3]>150&width(Peaks)>1000&width(Peaks)<5000)
#' Peaks <- Peaks[idx[1:500]]
#' }
#'
NULL

#' mm9-Genes
#'
#' Subset of Genes from the mm9 annotation that overlap with example Peaks in
#' the Cfp1-Peaks file.
#'
#' @name mm9-Genes
#' @docType data
#' @keywords data
#' @format contains, GR a GRanges object with 800 ranges
#' @usage data('mm9-Genes')
#' @examples
#' # data was created as follows:
#' \dontrun{
#' data('Cfp1-Peaks')
#' library(TxDb.Mmusculus.UCSC.mm9.knownGene)
#' txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene #shorthand (for convenience) txdb
#' GR <- transcripts(txdb)
#' ov <- findOverlaps(GR,Peaks)
#' GR <- GR[queryHits(ov)]
#' save(file = 'data/mm9-Genes.rData',GR)
#' }
#'
#'
NULL

#' DBAmmd Object for Cfp1 example
#' @usage data('MMD')
#' @name MMD
#' @docType data
#' @keywords data
#' @examples
#' # data was created as follows:
#' \dontrun{
#' library('MMDiff2')
#' library('MMDiffBamSubset')
#' # create metaData:
#' ExperimentData <- list(genome='BSgenome.Mmusculus.UCSC.mm9',
#'                       dataDir=system.file("extdata", package="MMDiffBamSubset"),
#'                       sampleSheet="Cfp1.csv")
#' MetaData <- list('ExpData' = ExperimentData)
#' MMD <- DBAmmd(MetaData)
#' data("Cfp1-Peaks")
#'
#' MMD <- setRegions(MMD,Peaks)
#' MMD <- getPeakReads(MMD,pairedEnd=FALSE, run.parallel=FALSE)
#'
#' MMD <- DBAmmd(MetaData)
#' MMD <- setRegions(MMD,Peaks)
#' MMD <- getPeakReads(MMD,pairedEnd=FALSE, run.parallel=FALSE)
#' MMD <- estimateFragmentCenters(MMD, shift=NULL, draw.on=FALSE)
#' MMD <- compHists(MMD, bin.length=20)
#' MMD <- compDists(MMD, dist.method = "MMD", run.parallel = FALSE)

#' group1 <- Samples(MMD)$Condition==1
#' names(group1) <- Samples(MMD)$SampleID
#' group2 <- Samples(MMD)$Condition==2
#' names(group2) <-  Samples(MMD)$SampleID
#' con <- list(group1=group1,
#'            group2=group2,
#'            name1='WT-Resc',
#'            name2='KO')
#'
#' MMD <- compPvals(MMD, contrasts=list(con))
#' }

NULL
