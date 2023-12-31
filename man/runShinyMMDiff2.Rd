% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runShinyMMDiff2.R
\name{runShinyMMDiff2}
\alias{runShinyMMDiff2}
\title{Shiny Application for interactive visualization of MMD,GMD and
Pearson Difference as well as plotting peaks}
\usage{
runShinyMMDiff2(MD, whichContrast = 1)
}
\arguments{
\item{MD}{DBAmmd Object. This Object can be created using \code{DBAmmd()}.}

\item{whichContrast}{index determining which of the set contrast should be used.
(DEFAULT: 1)}
}
\description{
Shiny Application for interactive visualization of MMD,GMD and
Pearson Difference as well as plotting peaks
}
\examples{
 if(interactive()){
  data("MMD")
runShinyMMDiff2(MMD)
}

}

