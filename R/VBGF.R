#' @title Von Bertalanffy Growth function (VBGF)
#
#' @description  This function applies the von Bertalanffy growth function. It
#'    allows to calculate ages from lengths based on the special, generalised
#'    or seasonalised VBGF.
#'
#' @param Linf infinite length for investigated species in cm, or
#' @param Winf infinite weight for investigated species in gramm
#' @param K growth coefficent for investigated species per year
#' @param t age
#' @param t0 theoretical time zero, at which individuals of this species hatch
#' @param b exponent of weight length relationship
#' @param D surface factor
#' @param L0 length at hatching for VBGF with L0
#' @param ts onset of the first oscillation relative to t0
#' @param C intensity of (sinusoid) growth oscillations
#'
#' @keywords function growth VBGF
#'
#' @examples
#' # with t0
#' t <- seq(0,6,0.1)
#' Lt <- VBGF(Linf=80, K=0.6, t=t, t0=-0.1)
#' plot(t, Lt, t="l")
#'
#' # with L0
#' t <- seq(0,6,0.1)
#' Lt <- VBGF(Linf=80, K=0.6, t=t, L0=2)
#' plot(t, Lt, t="l")
#'
#' #seasonalised
#' t <- seq(0,6,0.1)
#' Lt <- VBGF(Linf=80, K=0.6, t=t, t0=-0.1, ts=0.5, C=0.75)
#' plot(t, Lt, t="l")
#'
#' @details
#' Based upon which input parameters are given one of the following
#' VBGF types is applied: "special", "generalised", or "seasonalised" VBGF.
#'
#' @return A vector with estimated lengths corresponding to provided ages.
#'
#' @references
#'
#' Somers, I. F. (1988). On a seasonally oscillating growth function.
#' Fishbyte, 6(1), 8-11
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

VBGF <- function(t, Linf = NA, Winf = NA, K, t0 = 0,
                 b = 3, D = 1,
                 L0 = NA, ts = 0, C = 0){

  if(is.na(Linf) & is.na(Winf)) stop("You have to provide either Linf or Winf.")

  if(is.na(L0)){
    # generalised seasonalised VBGF for length
    if(is.na(Winf)) res <- Linf * (1 - exp(-K * D * (t - t0) -
                                            ((C*K*D)/(2*pi)) *
                                            ((sin(2*pi*(t-ts))) +
                                               (sin(2*pi*(t0-ts))))))^ (1/D)
    # generalised VBGF for weight
    if(is.na(Linf)) res <- Winf * (1 - exp(-K * D * (t - t0)))^(b/D)
  }

  # special VBGF for length with L0
  if(!is.na(L0)) res <- Linf - (Linf - L0) * exp(-K * t)

  return(res)
}

