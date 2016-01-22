#' @title Von Bertalanffy growth equation
#
#' @description  This is the von Bertalanffy growth equation with t0 or L0 respectively
#'
#' @param Linf infinite length for investigated species in cm [cm]
#' @param K growth coefficent for investigated species per year [1/year]
#' @param t time
#' @param t0 theoretical time zero, at which individuals of this species hatch
#'
#' @keywords function growth
#'
#' @examples
#' # with t0
#' t <- seq(0,6,0.1)
#' L <- growth_VB(Linf=80, K=0.6, t=t, t0=-0.1)
#' plot(t, L, t="l")
#'
#' with L0
#' t <- seq(0,6,0.1)
#' L <- growth_VB(Linf=80, K=0.6, t=t, L0=2)
#' plot(t, L, t="l")
#'
#'
#' @details
#'
#' @return A vector with lengths
#'
#' @references
#'
#' @export

growth_VB <- function(Linf, K, t, t0 = 0, L0 = NA){
  if(is.na(L0)) return(Linf * (1 - exp(-K * (t - t0))))
  if(!is.na(L0)) return(Linf - (Linf - L0) * exp(-K * t))
}


