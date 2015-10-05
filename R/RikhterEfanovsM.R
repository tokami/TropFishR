#' @title Rikhter and Efanov's formula
#
#' @description Rikhter and Efanov's formula to calculate natural mortality (M)
#'
#' @param Tm50 Age when 50% of the population is mature [year] ("the age of massive maturation")
#'
#' @examples
#' RikhterEfanovsM(Tm50 = 1)
#'
#' @details Rikhter and Efanov's formula to calculate natural mortality (M).
#'
#' @export

RikhterEfanovsM <- function(Tm50){

  M <- 1.521 / ( Tm50 ^ 0.720) - 0.155

  return(M)
}


