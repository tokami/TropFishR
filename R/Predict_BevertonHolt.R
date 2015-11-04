#' @title Beverton and Holt's prediction model
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param classes Midpoints of the length class as vector (length frequency data) or ages as vector (age composition data).
#' @param catch Catch as vector, or a matrix with catches of subsequent years if the catch curve with constat time intervals should be applied.
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch (default: 0).
#' @param Winf Infinite weight in grams!
#'
#' @examples
#'
#'
#' @details For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references
#'
#'
#' @export


# Nemipterus marginatus
K = 0.37
M = 1.1
Tc = 1.0
Tr = 0.4
t0 = -0.2
Winf = 286
#e.g.
F_PBH = 0.5

#Leiognathus spendens (Pauly 1980)
Wgama = 64
K = 1
t0 = -0.2
Tr = 0.2
M = 1.8

Tc = 0.2
Tc = 0.3
Tc = 1.0

############

#Xiphias gladius (Berkeley and Houde 1980)
Lgama = 309
K = 0.0949
M = 0.18

Lc = 118
Lc = 150



Predict_BevertonHolt <- function(Winf = NA, Linf = NA, K, t0 = NA,
                                  M, Tr = NA, Tc = NA, Lc = NA, FM = NA){

  if(is.na(Winf) & is.na(Linf)) stop("You have to provide Linf or Winf!")

  S.PBH <- exp(-K * (Tc - t0))
  Y_R.PBH <- F_PBH * (exp(-M*(Tc-Tr)) * Winf) * ((1/(F_PBH+M)) -
                                                    ((3*S.PBH)/(F_PBH + (M+K))) +
                                                    ((3*(S.PBH^2))/(F_PBH + (M+2*K))) -
                                                    ((S.PBH^3)/(F_PBH + (M+3*K))))

}


