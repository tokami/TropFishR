#' @title Heincke's method
#'
#' @description Method to estimate the total mortality (Z) from CPUE data.
#'
#' @param midLengths Midpoints of the length class as vector
#' @param catch Catch per sampling time as matrix or the total catch as vector.
#' @param Linf Infinite length for investigated species in cm [cm].
#'
#' @examples
#' data("ex.HeinckesZ")
#' HeinckesZ(cohort <- ex.HeinckesZ$cohort,
#'  age <- ex.HeinckesZ$age.t1,
#'  CPUE <- ex.HeinckesZ$CPUE)
#'
#' @details Z from CPUE data
#'
#' @references
#' Jones ???  Sparre?
#'
#' @export



HeinckesZ <- function(cohort, age, CPUE){

  df.HZ <- data.frame(cohort = cohort[1:(length(cohort)-1)],
                      age.t1 = age[1:(length(age)-1)],
                      CPUE = CPUE[1:(length(CPUE)-1)])

  result_Z <- list()
  for(i in 2:length(age)){
    Zi <- rep(NA,(length(age)-1))
    for(k in 1:(i-1)){
      Zi[k] <- round((1 / (age[i] - age[k])) *
                    (log(CPUE[k] /CPUE[i])),digits = 2)
    }
    result_Z[[i-1]] <- Zi
  }

  for(i in 1:length(result_Z)){
    df.HZ[[paste("Z",cohort[i+1],sep='_')]] <- result_Z[[i]]
  }

  print(df.HZ)
}
