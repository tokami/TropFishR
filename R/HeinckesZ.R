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
#'
#' @details Z from CPUE data
#'
#' @references
#' Jones ???  Sparre?
#'
#' @export


data("ex.HeinckesZ")

#HeinckesZ(
  cohort <- ex.CPUE_Z$cohort
  age <- ex.CPUE_Z$age.t1
  CPUE <- ex.CPUE_Z$CPUE

HeinckesZ <- function(cohort, age, CPUE){

  df.HZ <- data.frame(cohort = cohort,
                      age.t1 = age,
                      CPUE = CPUE)

  result_Z <- list()
  for(i in 2:length(age)){
    Zi <- round((1 / (age[i] - age[i-1])) *
                  (log(CPUE[i-1] /CPUE[i])),digits = 2)

  }


  for(i in 1:(length(CPUE)-1)){
    Zi <- round((1 / (age[i+1] - age[i])) * (log(CPUE[i] /CPUE[i+1])),digits = 2)
    delta_ti <- paste("t",i+1,"-t",i,sep='')
    result_Z[[i]] <- c(delta_ti,Zi)
  }

  result_Z <- do.call(rbind,result_Z)
  names(result_Z) <- c("time interval","Z")
  print(result_Z)

}
