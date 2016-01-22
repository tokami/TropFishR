#' @title Seasonalised version of the von Bertalanffy growth equation
#
#' @description  This is the seasonalised version of the von Bertalanffy growth equation
#'
#' @param Linf infinite length for investigated species in cm [cm]
#' @param K growth coefficent for investigated species per year [1/year]
#' @param t time
#' @param t0 theoretical time zero, at which individuals of this species hatch
#' @param ts summer point
#' @param C strength of oscillation
#'
#' @keywords function growth
#'
#' @examples
#' t <- seq(0,6,0.1)
#' L <- growth_soVB(Linf=80, K=0.6, t=t, t0=-0.1, ts=0.5, C=0.75)
#' plot(t, L, t="l")
#'
#'
#' @details
#'
#' @return A vector with lengths
#'
#' @references
#' Somers, I. F. (1988). On a seasonally oscillating growth function. Fishbyte, 6(1), 8-11
#'
#' @export

growth_soVB <- function(Linf, K, t, t0, ts, C){
  Linf * (1 - exp(-K * (t - t0) +
                    (((C * K)/(2 * pi)) *
                       sin(2 * pi * (t - ts))) - (((C * K)/(2 * pi)) *
                                                    sin(2 * pi * (t0 - ts)))))
}
