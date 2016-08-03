#' @name goatfish
#'
#' @title Yellowstriped goatfish data
#'
#' @description Data of Yellowstriped goatfish (\emph{Upeneus vittatus}) from Manila Bay,
#'    Philippines. Can be used for the estimation of the instantaneous mortality
#'    rate (Z) by means of \code{\link{catchCurve}}.
#'
#' @docType data
#' @format A list consisting of:
#' \itemize{
#'   \item \code{midLengths}: mid points of length classes,
#'   \item \code{catch}: a vector with catches in numbers,
#'   \item \code{Linf}: infinite length in cm [cm],
#'   \item \code{K}: growth coefficent per year [1/year].
#' }
#' @source Ziegler, B., 1979. Growth and mortality rates of some fishes of
#'   Manila Bay, Philippines, as estimated from analysis of length-frequencies.
#'   Thesis. Kiel University, 115 p.
#' @usage data(goatfish)
#' @keywords data dataset catchCurve length-frequency
#' @examples
#' data(goatfish)
#' str(goatfish)
#' summary(goatfish)
#'
#'
NULL
