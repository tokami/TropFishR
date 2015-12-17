#' @title Pauly's empirical formula for natural mortality (M)
#
#' @description Pauly's empirical formula to calculate natural mortality (M).
#'
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficient for investigated species per year [1/year].
#' @param Temp Average annual temperature at the surface in degrees centigrade.
#' @param schooling Logical parameter accounting for schooling behaviour, according to Pauly (1983). Default is FALSE.
#'
#' @examples
#' M_Pauly(Linf = 80, K = 0.5, temp = 25)
#'
#' @details Pauly's empirical formula to calculate natural mortality (M). If accounting for schooling behaviour M is multiplied by 0.8 according to Pauly (1983).
#'
#' @references
#' Pauly, D., 1980. On the interrelationships between natural mortality, growth parameters,
#' and mean environmental temperature in 175 fish stocks. \emph{J.Cons.CIEM}, 39(2):175-192
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

M_Pauly <- function(Linf, K, temp, schooling = FALSE){

  M <- exp( -0.0152 - 0.279 * log(Linf) + 0.6543 * log(K) + 0.463 * log(temp))

  if(schooling == TRUE){
    M <- 0.8 * M
  }

  return(M)
}


