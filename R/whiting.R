#' @name whiting
#'
#' @title Whiting data
#'
#' @description A dataset of North Sea whiting \emph{Merlangius merlangus}
#'    caught during the period 1974-1980. Can be used for \code{\link{catchCurve}}
#'    and \code{\link{VPA}}.
#'
#' @docType data
#' @format A list consisting of:
#' \itemize{
#'   \item \code{age}: a vector with age groups,
#'   \item \code{M}: natural mortality rate,
#'   \item \code{a}: length-weight relationship coefficent (W = a * L^b),
#'   \item \code{b}: length-weight relationship coefficent (W = a * L^b),
#'   \item \code{catch}:a matrix with catches from 1974 to 1980.
#' }
#' @source ICES, 1981. Report of the \emph{Ad hoc} working group on the use of effort data
#'    in assessment, Copenhagen, 2-6 March 1981. \emph{ICES C.M.} 1981/G:5 (mimeo)
#' @usage data(whiting)
#' @keywords dataset data VPA catchCurve
#' @examples
#' data(whiting)
#' str(whiting)
#' summary(whiting)
#'
#'
NULL
