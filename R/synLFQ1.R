#' @name synLFQ1
#'
#' @title Synthetic length-frequency data I
#'
#'
#' @description Synthetic length-frequency data as provided in Sparre & Venema (1998).
#'    Can be used to apply the function \code{\link[TropFishR]{Bhattacharya}}.
#'
#' @docType data
#'
#' @format A list consisting of:
#' \itemize{
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column).
#' }
#'
#' @source Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#'    Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @usage data(synLFQ1)
#' @keywords data dataset length-frequency
#'
#' @examples
#' data(synLFQ1)
#' str(synLFQ1)
#' summary(synLFQ1)
#'
#'
NULL
