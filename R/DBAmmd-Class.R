#' Class DBAmmd
#'
#' Class \code{DBAmmd} defines a container for differential binding analysis
#' using MMDiff2. For this class a number of methods is foreseen, among which
#' accessors for every slot.
#' As MetaData it needs to contain the path to the data directory and the
#' name of a sampleSheet csv file. Additionally, the Regions of interests must
#' be provided as GRanges Object.
#'
#' @section Slots:
#' \describe{
#' \item{\code{MetaData}:}{MetaData}
#' \item{\code{rowRanges}:}{GRanges list object containing Regions of Interests
#' (Peaks)}
#' \item{\code{Reads}:}{List containing positions of mapped reads, i.e. exact
#' start and end positions of mapped fragments, not coverage. In the case of
#' single-end reads, the left most postions of fragments mapping to the positive
#' strands and the right most positions of fragments
#' mapping to the negative strands are stored in "Left.p" and "Right.n".
#' Use \code{getPeakReads} to fill this slot and \code{estimateFragmentCenters}
#' to add the (estimated) positions of fragment centers.}
#' \item{\code{RawTotalCounts}:}{ m x n matrix containing total counts of reads
#' mapping
#' to m peaks and in n samples (including input samples)}
#' \item{\code{RawCounts.p}:}{ m x n matrix containing counts of reads mapping
#' to positive (forward) strand}
#' \item{\code{RawCounts.n}:}{ m x n matrix containing counts of reads mapping
#' to negative (forward) strand}
#' \item{\code{Hists}:}{ For each Peak matrix}
#' \item{\code{DISTs}:}{ for each Peak computed distances between pairs of
#' samples}
#' \item{\code{mCounts}:}{ (for internal use only, counts averaged over pairs of
#' samples)}
#' \item{\code{Contrasts}:}{ 'list'}
#' }
#'
#' @note For this class,
#' @seealso \code{\link{DBAmmd-Accessors}}
#' @return DBAmmd Object
#'
#' @importFrom utils read.csv
#' @name DBAmmd-class
#' @rdname DBAmmd-class
#' @exportClass DBAmmd
#' @aliases DBAmmd
#' @author Gabriele Schweikert
#'


setClass(Class="DBAmmd", slots=c(
  MetaData="list",
  rowRanges="GRanges",
  Reads='list',   #contains lists of read start positions
  RawTotalCounts='matrix',
  RawCounts.p='matrix',
  RawCounts.n='matrix',
  Hists='list',
  DISTs='list',
  mCounts='matrix',
  Contrasts='list')
)

## Constructor method

#'@export
DBAmmd <- function(MetaData=NULL){
  rowRanges <- GRanges()
  Reads <- list()
  RawTotalCounts <- matrix(nrow = 0, ncol = 0)
  RawCounts.p <- matrix(nrow = 0, ncol = 0)
  RawCounts.n <- matrix(nrow = 0, ncol = 0)
  Hists <- list()
  DISTs <- list()
  mCounts <- matrix(nrow = 0, ncol = 0)
  Contrasts <- list()

  if (is.null(MetaData)){
    MetaData=list()}
  else {
    if (is.null(MetaData$ExpData))
      stop('MetaData requires field ExpData
           (list containing dataDir and SampleSheet).')
    if (is.null(MetaData$ExpData$dataDir) |
        is.null(MetaData$ExpData$sampleSheet))
      stop('MetaData requires field ExpData
           (list containing dataDir and SampleSheet).')
    oldwd <- setwd(MetaData$ExpData$dataDir)
    samples <- read.csv(MetaData$ExpData$sampleSheet,
                        header = TRUE,colClasses = "character")
    setwd(oldwd)
    MetaData$ExpData$samples=samples
  }

  new("DBAmmd",
      MetaData=MetaData,
      rowRanges=rowRanges,
      Reads=Reads,
      RawTotalCounts=RawTotalCounts,
      RawCounts.p=RawCounts.p,
      RawCounts.n=RawCounts.n,
      Hists=Hists,
      DISTs=DISTs,
      mCounts=mCounts,
      Contrasts=Contrasts)
}
