#' @title Methods to calculate Z from CPUE data
#'
#' @description Method to estimate the total mortality (Z) from CPUE data.
#'
#' @param cohort bla bla
#' @param age bla bla
#' @param CPUE bla bla
#' @param method One of 3 methods to caluclate Z: "standard", "Heincke"
#'
#' @examples
#' data("ex.Z_CPUE")
#' with(ex.Z_CPUE,Z_CPUE(cohort, age.t1, CPUE, method = "standard"))
#'
#' @details Z from CPUE data
#'
#' @return A matrix with total moratlity rates with column names ... and row names ... .
#'
#' @references
#' Jones ???  Sparre?
#'
#' @export

Z_CPUE <- function(cohort, age, CPUE, method){

  if(method == "standard"){
    df.HZ <- data.frame(cohort = cohort[1:(length(cohort)-1)])

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
      df.HZ[[paste(cohort[i+1])]] <- result_Z[[i]]
    }
    print(df.HZ)
  } else if(method == "Heincke"){
    CPUE[4] <- sum(CPUE[4:length(CPUE)])
    CPUE.H.n <- CPUE[2:length(CPUE)]
    CPUE.H.d <- CPUE[1:length(CPUE)]
    Z.H = - log( (sum(CPUE.H.n)) /
                 (sum(CPUE.H.d)))
    print(paste("Total mortality according to Heincke's method: ",Z.H))
  }
}
