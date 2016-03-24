#' @title Pauly's empirical formula for the estimation of M
#
#' @description Function to calculate the instantaneous natural mortality rate (M)
#'      according to Pauly's empirical formula.
#'
#' @param Linf infinite length for investigated species in cm [cm].
#' @param K growth coefficient for investigated species per year [1/year].
#' @param temp average annual temperature at the surface in degrees centigrade.
#' @param schooling logical; if TRUE it is accounted for the schooling behaviour of
#'      the species. Default is FALSE.
#'
#' @keywords function mortality M
#'
#' @examples
#' M_Pauly(Linf = 80, K = 0.5, temp = 25)
#'
#' @details If accounting for schooling behaviour M is multiplied by 0.8 according to
#'    Pauly (1983).
#'
#' @return A numeric value representing M.
#'
#' @references
#' Pauly, D., 1980. On the interrelationships between natural mortality, growth parameters,
#' and mean environmental temperature in 175 fish stocks. \emph{J.Cons.CIEM}, 39(2):175-192
#'
#' Pauly, D., 1983. Some simple methods for the assessment of tropical fish stocks.
#' \emph{FAO Fish.Tech.Pap.}, (234): 52p. Issued also in French and Spanish
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

M_Pauly <- function(Linf, K, temp, schooling = FALSE){

  M <- exp( -0.0152 - 0.279 * log(Linf) + 0.6543 * log(K) + 0.463 * log(temp))

  if(schooling == TRUE){
    M <- 0.8 * M
  }

  return(M)
}


