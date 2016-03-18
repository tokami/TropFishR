#' @title Rikhter and Efanov's formula for the estimation of M
#'
#' @description Function to estimate the instantaneous natural mortality rate (M)
#'      according to Rikhter and Efanov's formula.
#'
#' @param tm50 age when 50\% of the population is mature [year]
#'      ("age of massive maturation").
#'
#' @keywords function mortality M
#'
#' @examples
#' M_RikhterEfanov(tm50 = 1)
#'
#' @return A numeric value representing M.
#'
#' @references
#' Rikhter, V.A., and V.N. Efanov, 1976. On one of the approaches to estimation of natural
#' mortality of fish populations. \emph{ICNAF Res.Doc.}, 76/VI/8: 12p.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

M_RikhterEfanov <- function(tm50){

  M <- 1.521 / ( tm50 ^ 0.720) - 0.155

  return(M)
}


