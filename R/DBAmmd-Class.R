#' Class DBAmmd
#'
#' The \code{DBAmmd} Class defines a container for differential binding analysis
#' using MMDiff2. For this class a number of methods is foreseen, among which
#' accessors for every slot.
#' As MetaData, it needs to contain the path to the data directory and the
#' name of a sampleSheet csv file.
#'
#' @section Constructor:
#' \code{DBAmmd()}{returns an empty DBAmmd Object.} \cr
#' \code{DBAmmd(MetaData)}{ initializes a DBAmmd Object for a new
#' Experiment. \cr
#' (See below and the package vignette for more details.)}
#'
#' @section Slots:
#' \describe{
#' \item{\code{MetaData}:}{List containing an \code{ExpData} and an
#' \code{AnaData} compartment. "ExpData" needs a \code{dataDir} and a
#' \code{SampleSheet} entry. A \code{genome} entry, which should be a valid
#' \code{BSGenome} name, is useful to find sequence motifs. (Note the genome
#' version needs to correspond to the one used for the read alignment.
#' Use \code{available.genomes()} to find the right name.) The \code{AnaData}
#' entry is used to store and access parameters for the MMDiff2 Analysis, like the sigma
#' of the RBF Kernel.}
#' \item{\code{rowRanges}:}{GRanges object containing Regions of Interests
#' (Peaks)}
#' \item{\code{Reads}:}{List containing positions of mapped reads, i.e. exact
#' start and end positions of mapped fragments. In the case of
#' single-end reads, the left most postions of fragments mapping to the positive
#' strands and the right most positions of fragments
#' mapping to the negative strands are stored in "Left.p" and "Right.n".
#' Use \code{getPeakReads} to fill this slot and \code{estimateFragmentCenters}
#' to add the (estimated) positions of fragment centers.}
#' \item{\code{RawTotalCounts}:}{ m x n matrix containing total counts of reads
#' mapping
#' to m peaks in n samples (including input samples)}
#' \item{\code{RawCounts.p}:}{ m x n matrix containing counts of reads mapping
#' to positive (forward) strand}
#' \item{\code{RawCounts.n}:}{ m x n matrix containing counts of reads mapping
#' to negative (reverse) strand}
#' \item{\code{Hists}:}{ List of lists, each of length m (number of Peaks).
#' Compartments could be 'Left.p','Right.n','Left.n','Right.p','Center.n',
#' 'Center.p','Center','Left','Right', defining whether left or right ends or
#' centers of fragments should be considered for positive ('p') or negative ('n')
#' strand, or both strands combined. For a given compartment there is one entry per
#' peak, which is a n x L_i matrix, where n is the number of samples and L_i is
#' the number of bins used to cover the extend of the peak. Note, L_i varies
#' between peaks of different lengths. See \code{compHists()} for more details.}
#' \item{\code{DISTs}:}{List with compartments for different methods to compute
#' distances (e.g. MMD). Each compartment contains a m x N matrix with computed
#' distances for each Peak between N pairs of samples.
#' See \code{compDists()} for more details.}
#' \item{\code{mCounts}:}{ (for internal use only)}
#' \item{\code{Contrasts}:}{List of lists.
#' Each entry contains a contrast i.e. the definition of two groups that should
#' be compared to each other in a differential analysis. A Contrast needs entries
#'  "name1", "name2" for group names, as well as group memberships given in
#' "group1" and "group2". Results of a differential test for this contrast are
#' stored in an entry given by the method name, e.g. "MMD.locfit"}
#' }
#'
#'
#' @importFrom utils read.csv
#'
#' @seealso \code{\link{DBAmmd-Accessors}},\code{\link{getPeakReads}}
#'
#'
#' @return DBAmmd Object
#'
#' @name DBAmmd-class
#' @rdname DBAmmd-class
#' @exportClass DBAmmd
#' @aliases DBAmmd
#' @author Gabriele Schweikert
#'
#' @examples
#'
#' ## Example using a small data set provided in the MMDiffBamSubset package
#'
#' # setting the Experiment meta data:
#' ExpData <- list(dataDir=system.file("extdata", package="MMDiffBamSubset"),
#'            sampleSheet="Cfp1.csv")
#'
#' MetaData <- list('ExpData' = ExpData)
#'
#' # Creating a DBAmmd data set:
#' MMD <- DBAmmd(MetaData)
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
