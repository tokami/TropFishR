% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.Bhattacharya.R
\name{plot.Bhattacharya}
\alias{plot.Bhattacharya}
\title{Bhattacharya plot}
\usage{
\method{plot}{Bhattacharya}(x, analysisPlot = TRUE, ...)
}
\arguments{
\item{x}{a list of the class \code{"Bhattacharya"} containing the results of
Bhattacharya's method.}

\item{analysisPlot}{logical; indicating wheter the anaylsis graph with the regression
lines should be created}

\item{...}{additional options of the \code{\link{plot}} function}
}
\description{
This function plots the seperated frequency distributions
     and selected regression lines of \link[TropFishR]{Bhattacharya} method.
}
\details{
This function plots the results of the Bhattacharya method.
}
\examples{
\dontrun{
 data(synLFQ1)
 output <- Bhattacharya(param = synLFQ1)
 plot(output)
}

}
\references{
Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
}
