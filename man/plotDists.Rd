% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotDists.R
\name{plotDists}
\alias{plotDists}
\title{plotDists}
\usage{
plotDists(MD, dist.method = "MMD", whichContrast = 1, which.group1 = NULL,
  which.group2 = NULL, diff.method = "MMD.locfit", bUsePval = FALSE,
  th = 0.1, title = NULL, what = 3, xlim = NULL, ylim = NULL,
  xlog10 = TRUE, Peak.IDs = NULL, withLegend = TRUE,
  shiny_df_opt = FALSE)
}
\arguments{
\item{MD}{DBAmmd Object. This Object can be created using \code{DBAmmd()}.}

\item{dist.method}{specify method used for distances
between samples. Currently only Maximum Mean Discrepancy (MMD)
and Kolmogorov-Smirnov (KS) implemented.
(DEFAULT: 'MMD')}

\item{whichContrast}{index determining which of the set contrast should be used.
(DEFAULT: 1)}

\item{which.group1}{subset samples from group1 (DEFAULT: NULL)}

\item{which.group2}{subset samples from group2 (DEFAULT: NULL)}

\item{diff.method}{which method to use to determine significant peaks
(DEFAULT: 'MMD.locfit')}

\item{bUsePval}{if TRUE p-values instead of FDRs are used (DEFAULT: FALSE)}

\item{th}{significance threshold for differential called peaks (DEFAULT: 0.1)}

\item{title}{an overall title for the plot (DEFAULT: NULL)}

\item{what}{which dists to overlay: 1: only between group distances,
2: between and within group distances, 3: between and within group distances,
and significant peaks highlightend   (DEFAULT: 3)}

\item{xlim}{specify x range (DEFAULT: NULL)}

\item{ylim}{specify y range (DEFAULT: NULL)}

\item{xlog10}{should x range be plotted in log10 scale (DEFAULT: TRUE)}

\item{Peak.IDs}{Highlight specific subset of peaks (DEFAULT: NULL)}

\item{withLegend}{(DEFAULT: TRUE)}

\item{shiny_df_opt}{Option returns a dataframe for shiny (DEFAULT: FALSE)}
}
\description{
scatterplot showing distances between peaks
}
\examples{
data("MMD")
plotDists(MMD, whichContrast=1)

}

