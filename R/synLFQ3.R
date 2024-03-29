#' @name synLFQ3
#'
#' @title Synthetic length frequency data III
#'
#' @description Synthetic length frequency data from Sparre & Venema (1998). Can be
#'    used for the estimation of the instantaneous total mortality rate (Z) by means
#'    of \code{\link{Z_BevertonHolt}}.
#'
#' @docType data
#'
#' @format A list consisting of:
#' \itemize{
#'   \item \code{midLengths}: midpoints of the length classes,
#'   \item \code{Linf}: infinite length for investigated species in cm [cm],
#'   \item \code{K}: growth coefficent for investigated species per year [1/year],
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \code{catch}: a vector with catches.
#' }
#'
#' @source Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
#'      assessment. Part 1. Manual. \emph{FAO Fisheries Technical Paper},
#'      (306.1, Rev. 2). 407 p.
#'
#' @usage data(synLFQ3)
#' @keywords data dataset length-frequency
#'
#' @examples
#' data(synLFQ3)
#' str(synLFQ3)
#' summary(synLFQ3)
#'
#'
NULL
