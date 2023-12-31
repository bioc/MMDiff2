% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DBAmmd-Accessors.R
\docType{methods}
\name{DBAmmd-Accessors}
\alias{Contrast,DBAmmd-method}
\alias{Counts,DBAmmd-method}
\alias{DBAmmd-Accessors}
\alias{Dists,DBAmmd-method}
\alias{Genome,DBAmmd-method}
\alias{Hists,DBAmmd-method}
\alias{Reads,DBAmmd-method}
\alias{Regions,DBAmmd-method}
\alias{Samples,DBAmmd-method}
\alias{metaData,DBAmmd-method}
\alias{numPeaks,DBAmmd-method}
\alias{numSamples,DBAmmd-method}
\alias{setContrast,DBAmmd-method}
\alias{setRegions,DBAmmd-method}
\title{Extract data from DBAmmd objects}
\usage{
\S4method{Genome}{DBAmmd}(x)

\S4method{Samples}{DBAmmd}(x)

\S4method{numPeaks}{DBAmmd}(x)

\S4method{numSamples}{DBAmmd}(x)

\S4method{metaData}{DBAmmd}(x)

\S4method{Regions}{DBAmmd}(x)

\S4method{Reads}{DBAmmd}(x, whichPos = "Center")

\S4method{Counts}{DBAmmd}(x, whichCounts = "T")

\S4method{Hists}{DBAmmd}(x, whichPos = "Center")

\S4method{Dists}{DBAmmd}(x, dist.method = NULL)

\S4method{Contrast}{DBAmmd}(x, whichContrast = 1)

\S4method{setRegions}{DBAmmd}(x, Regions)

\S4method{setContrast}{DBAmmd}(x, contrast)
}
\arguments{
\item{x}{a DBAmmd Object. An empty instance can be created using \code{DBAmmd()}.
(See \code{\link{DBAmmd-class}} for more details.)}

\item{whichPos}{specifies which relative positions of mapped fragments
should to be considered.
Can be one of: 'Left.p', 'Right.p', 'Right.p' and 'Left.n':
Start and end positions of fragments mapping to positive or negative strand,
respectively ('Right.p' and 'Left.n' are not available for single-end reads).
Additionally inferred positions: 'Center.n','Center.p','Center','Left','Right'.
(DEFAULT: 'Center')}

\item{whichCounts}{can be 'T': total counts, or
'p','n': counts of reads mapping to positive, negative strand, respectively.}

\item{dist.method}{specify method used for distances
between samples. Currently only Maximum Mean Discrepancy (MMD)
and Kolmogorov-Smirnov (KS) implemented.
(DEFAULT: 'MMD')}

\item{whichContrast}{index determining which of the set contrast should be used.
(DEFAULT: 1)}

\item{Regions}{GRanges Object specifying the Regions of Interesst / Peaks.}

\item{contrast}{determines how to set a new contrast for differential analysis.
A contrast can be automatically created either 'byCondition', or 'byTissue'.
The Contrast can also be manually set (see vignette for details).}
}
\value{
\code{Genome(x)} returns the name of the used genome version, if set
in the metaData.

\code{Samples(x)} returns the information which was provided in the
SampleSheet.csv to describe the data.

\code{numPeaks(x)} returns the number of Peaks / Regions of Interest
that are associated with the DBAmmd object.

\code{numSamples(x)} returns the number of samples associated with the
DBAmmd object.

\code{metaData(x)} returns the metaData associated with the
DBAmmd object.

\code{Regions(x)} returns the Peaks / Regions of Interest that are
associated with the DBAmmd object.

\code{Reads(x,whichPos)} returns the Reads mapping to the Regions of Interest.

\code{Counts(x,whichCounts)} returns a m x n matrix containing the
Counts of Reads mapping to the Peaks / Regions of Interest.
Depending on the value of 'whichCounts', total counts ('T'),
or counts of reads mapping to positive ('p'), or negative strand ('n')
are returnt. See \code{\link{getPeakReads}} for more details.

\code{Hists(x,whichPos)} returns a list of matrices of length m
(number of Peaks). Each matrix is a n x L_i matrix, where n is the number of
samples and L_i is the number of bins used to cover
the extend of the peak. Note, L_i varies between peaks of different lengths.
See \code{\link{compHists}} for more details.

\code{Dists(x,dist.method)} returns a matrix containing distances
between pairs of samples for each peak. See \code{\link{compDists}} for
more details.

\code{Contrast(x,whichContrast)} returns the specified contrast.

\code{setRegions(x,Regions)} returns a DBAmmd Object with set
Peaks / Regions of Interests.

\code{setContrast(x,contrast)} returns a DBAmmd Object
with a set contrast.
}
\description{
This help file describes different ways to access the slots and values
contained in a \code{\link{DBAmmd-class}} objects.
}
\examples{
data("MMD")

Samples(MMD)
Genome(MMD)
numPeaks(MMD)
numSamples(MMD)
metaData(MMD)
R <- Regions(MMD)
Pos <- Reads(MMD)
C <- Counts(MMD)
H <- Hists(MMD)
D <- Dists(MMD)
C1 <- Contrast(MMD)

}
\seealso{
\code{\link{DBAmmd-class}}
}

