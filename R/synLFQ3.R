#' @name synLFQ3
#'
#' @title Synthetic length frequency data III from Sparre & Venema (1998)
#'
#'
#' @description Synthetic length frequency data III. Can be used for the calculation of Z
#'    by means of \code{\link{BevertonHoltsZ}}.
#'
#' @docType data
#'
#' @format A list consisting of:
#' \itemize{
#'   \item \code{$midLengths} a vector of the mid lengths of the length groups,
#'   \item \code{$Linf} Infinite length for investigated species in cm [cm],
#'   \item \code{$K} Growth coefficent for investigated species per year [1/year],
#'   \item \code{$t0} Theoretical time zero, at which individuals of this species hatch,
#'   \item \code{$catch} a vector with catches.
#' }
#'
#' @source Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment. Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @usage data(synLFQ3)
#' @keywords data dataset length-frequency
#'
#' @examples
#' data(synLFQ3)
#' str(synLFQ3)
#'
#'
NULL
