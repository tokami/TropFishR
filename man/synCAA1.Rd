% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/synCAA1.R
\docType{data}
\name{synCAA1}
\alias{synCAA1}
\title{Synthetic Catch-at-age data I}
\format{
A list consisting of:
\itemize{
  \item \code{midAge}: a vector of the mid ages of the age groups,
  \item \code{catch}: a matrix with the catches for different years.
}
}
\source{
Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
   Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
}
\usage{
data(synCAA1)
}
\description{
Synthetic Catch-at-age data from Sparre & Venema (1998). Can be used
   for the estimation of the instantaneous total mortality rate (Z) by means
   of \code{\link{Z_BevertonHolt}}.
}
\examples{
data(synCAA1)
str(synCAA1)
summary(synCAA1)


}
\keyword{CAA}
\keyword{Catch-at-age}
\keyword{age}
\keyword{data}
\keyword{dataset}
