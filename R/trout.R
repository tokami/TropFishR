#' @name trout
#' @title Trout data
#'
#' @description This dataset contains information about coral trout
#'  (\emph{Plectropomus leopardus}) caught near Heron Island (Great Barrier Reef,
#'  Australia) in October 1971.
#'  It can be used for the estimation of growth parameters, \code{\link{ELEFAN}}.
#'
#' @docType data
#'
#' @format A list consisting of:
#' \itemize{
#'   \item \strong{sample.no} sample numbers,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column).
#' }
#'
#' @source Goeden, G.B., 1978. A monograph of the coral trout \emph{Plectropomus leopardus}
#' (Lacepede). \emph{Res.Bull.Fish.Serv.Queensl.}, (1):42 p.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @usage data(trout)
#' @keywords data dataset trout length-frequency
#' @examples
#'
#' data(trout)
#' str(trout)
#' summary(trout)
#'
#'
NULL
