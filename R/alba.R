#' @name alba
#'
#' @title Length-frequency data of the clam Abra alba
#'
#' @description Length-frequency data of the clam Abra alba as presented
#' by Brey et al. (1988). Includes 7 approximately bi-monthly samplings of
#' A. alba binned into 14 length classes. The data is of class "lfq", which
#' can be used used e.g. in functions estimating growth parameters:
#' \code{\link[TropFishR]{ELEFAN}}, \code{\link[TropFishR]{ELEFAN_SA}} or
#' \code{\link[TropFishR]{ELEFAN_GA}}.
#'
#'
#' @docType data
#'
#' @format A list of class lfq consisting of:
#' \itemize{
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column).
#' }
#'
#' @source Brey, T., Soriano, M., and Pauly, D. 1988.
#' Electronic length frequency analysis: a revised and
#' expanded user's guide to ELEFAN 0, 1 and 2.
#' ICLARM Contrib. 261.
#'
#' @usage data(alba)
#' @keywords data dataset length-frequency
#'
#' @examples
#' data(alba)
#'
#' # plot raw catch frequencies
#' plot(alba, Fname = "catch")
#'
#' # plot restructured frequencies
#' alba <- lfqRestructure(alba, MA=5)
#' plot(alba, Fname = "rcounts")
#'
#' # growth parameters: see the examples of ELEFAN_SA and ELEFAN_GA
#'
NULL
