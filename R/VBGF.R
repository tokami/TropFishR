#' @title Von Bertalanffy Growth function (VBGF)
#
#' @description  This function applies the von Bertalanffy growth function (VBGF).
#'    It allows to calculate ages from lengths or lengths from ages based on the special,
#'    generalised or seasonalised VBGF.
#'
#' @param t ages for which to calculate corresponding lengths, or
#' @param L lengths for which to calculate corresponding ages
#' @param Linf infinite length for investigated species in cm, or
#' @param Winf infinite weight for investigated species in gramm
#' @param K growth coefficent for investigated species per year
#' @param t0 theoretical time zero, at which individuals of this species hatch (default: 0)
#' @param b exponent of weight length relationship (default: 3)
#' @param D surface factor (default: 1)
#' @param L0 length at hatching for VBGF with L0
#' @param ts onset of the first oscillation relative to t0
#' @param C intensity of (sinusoid) growth oscillations. Default is no oscillation (C = 0).
#'
#' @keywords function growth VBGF
#'
#' @examples
#' # calculation of lengths
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
#' # with Winf
#' t <- seq(0,6,0.1)
#' Wt <- VBGF(Winf=4000, K=0.8, t=t)
#' plot(t, Wt, t="l")
#'
#' # seasonalised VBGF
#' t <- seq(0,6,0.1)
#' Lt <- VBGF(Linf=80, K=0.6, t=t, t0=-0.1, ts=0.5, C=0.75)
#' plot(t, Lt, t="l")
#'
#'
#' # calculation of ages
#' L <- seq(2,200,0.1)
#' t <- VBGF(L = L, Linf=210, K=0.8, C= 0.75)
#' plot(t, L, t="l")
#'
#'
#' @details
#' Based upon which input parameters are given one of the following
#' VBGF types is applied: "special", "generalised", or "seasonalised" VBGF.
#'
#' @return A vector with estimated lengths corresponding to provided ages.
#'
#' @references
#' Somers, I. F. (1988). On a seasonally oscillating growth function.
#' Fishbyte, 6(1), 8-11
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

VBGF <- function(t = NA, L = NA, Linf = NA, Winf = NA, K, t0 = 0,
                 b = 3, D = 1,
                 L0 = NA, ts = 0, C = 0){

  if(is.na(Linf) & is.na(Winf)) stop("You have to provide either Linf or Winf.")
  if(is.na(L[1]) & is.na(t[1])) stop("Please provide a vector with either potential ages or lengths.")

  if(is.na(L0)){
    # generalised seasonalised VBGF for length
    if(is.na(Winf) & is.na(L[1])) res <- Linf * (1 - exp(-K * D * (t - t0) + (((C*K*D)/(2*pi)) * sin(2*pi*(t-ts))) - (((C*K*D)/(2*pi)) * sin(2*pi*(t0-ts))))) ^ (1/D)

    # OLD:
      # res <- Linf * (1 - exp(-K * D * (t - t0) -
      #                                       ((C*K*D)/(2*pi)) *
      #                                       ((sin(2*pi*(t-ts))) +
      #                                          (sin(2*pi*(t0-ts))))))^ (1/D)


    if(is.na(Winf) & is.na(t[1])){
      if(D == 1 & C == 0) res <- t0 - (log(1-L/Linf)/K)
      # lookup table for soVBGF
      if(D != 1 | C != 0){
        tmax <- (t0 * K * D - log(1 - exp(D*log(L[length(L)]/Linf)))) / (K*D)
        if(is.na(tmax)) tmax <- 40
        lookup_age <- seq(0,(tmax+10),0.001)
        lookup_length <-  Linf * (1 - exp(-K * D * (lookup_age - t0) + (((C*K*D)/(2*pi)) *
                                                                          sin(2*pi*(lookup_age-ts))) -
                                            (((C*K*D)/(2*pi)) * sin(2*pi*(t0-ts))))) ^ (1/D)

        # OLD:
        # lookup_length <-  Linf * (1 - exp(-K * D * (lookup_age - t0) -
        #                                     ((C*K*D)/(2*pi)) *
        #                                     ((sin(2*pi*(lookup_age-ts))) +
        #                                        (sin(2*pi*(t0-ts))))))^ (1/D)

        lookup_ind <- lapply(X = L, FUN = function(x) which.min((lookup_length - x)^2))
        res <- lookup_age[unlist(lookup_ind)]
      }
    }


    # generalised VBGF for weight
    if(is.na(Linf) & is.na(L[1])) res <- Winf * (1 - exp(-K * D * (t - t0)))^(b/D)
    if(is.na(Linf) & is.na(t[1])) res <- (D * K * t0 - log(1 - exp((D * log(L/Winf))/(b)))) / (K*D)
  }

  # special VBGF for length with L0
  if(!is.na(L0) & is.na(L[1])) res <- Linf - (Linf - L0) * exp(-K * t)
  if(!is.na(L0) & is.na(t[1])) res <- log((L - Linf) / -(Linf - L0)) / -K

  return(res)
}

