#' @title  Production models with time series fitting
#'
#' @description This function applies the production models under non-equilibrium
#'      conditions by applying time series fitting using non-linear least squares
#'      minimisation.
#'
#' @param data a dataframe of parameters
#' \itemize{
#'   \item \strong{year} years,
#'   \item \strong{yield} catch in weight of fishery per year,
#'   \item \strong{effort} fishing effort per year,
#'   \item \strong{CPUE} catch per unit of effort per year (optional).
#' }
#' @param method indicating if Schaefer or Fox model should be applied. First assumes a
#'    logistic relationship between growth rate and biomass, whereas second assumes it to
#'    foolow the Gompertz distribution (Richards 1959). Default is the dynamic Schaefer model.
#' @param SSE_method method which is used to calculate the error sum of squares of
#'    the observed and predicted CPUE values in the non-linear minimalisation approach.
#'    Options are "standard", for \code{sum((CPUE - predicted CPUE)^2)}, "log" for
#'    \code{sum((log(CPUE) - log(predicted CPUE))^2)}, and "Thiels_U_statistic" for
#'    \code{sqrt(sum(CPUE - predicted CPUE)/sum(CPUE(t) - CPUE(t-1)))}.
#' @param B0_init numerical; if realistic initial estimate for virgin biomass is available.
#'    If NA initial estimate for virgin biomass is set to two times average yield of all
#'    or part of yield values (see \code{B0_est}).
#' @param B0_est intital value of virgin biomass estimating using all yield values (NA) or
#'    first years of time series, then provide numerical representing number of years
#' @param effort_unit multiplication factor for the unit of effort. Default is 1.
#' @param plot logical; if TRUE (default) a graph is displayed
#'
#' @keywords function biomass MSY
#'
#' @examples
#' data(emperor)
#' #x2 <- prod_mod_ts(data = emperor, method = "Fox", SSE_method = "standard")
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{Bvec}: biomass vector,
#'   \item \strong{CPUE_hat}: predicted CPUE,
#'   \item \strong{K}: carrying capacity,
#'   \item \strong{r}: population growth rate,
#'   \item \strong{q}: catchability coefficient,
#'   \item \strong{MSY}: maximum sustainabale yield (MSY),
#'   \item \strong{Bmsy}: biomass yielding in MSY,
#'   \item \strong{Emsy}: fishing effort yielding in MSY
#'   \item \strong{Fmsy}: fishing mortality yielding in MSY,
#' }
#'
#' @details Either catch per unit of effort (CPUE) is inserted
#'    into the model directly (by a column \strong{CPUE}) or CPUE is calculated from
#'    the catch and effort, then these two vectors should have required units.
#'    Whenever a good estimate for the virigin biomass is available, this estimate
#'    should be inserted for \code{B_init}. The default approach for the initial
#'    estimate of the virgin biomass is to multiply the average yield by 2 (Dharmendra
#'    and Solmundsson, 2005). Alternatively, just a part of the time series of
#'    yield values can be choosen to represent the virgin biomass.
#'    With the Fox model it is recommeneded to take Thiels_U_statistic as the method
#'    to calculate SSE.
#'
#' @references
#' Dharmendra, D., Solmundsson, J., 2005. Stock assessment of the offshore Mauritian banks
#' using dynamic biomass models and analysis of length frequency of the Sky Emperor
#' (\emph{Lethrinus mahsena}). Fisheries Training Program The United Nations University, 61
#'
#' Hilborn, R. and Walters, C., 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York
#'
#' Richards, F. J. (1959). A flexible growth function for empirical use.
#' \emph{Journal of experimental Botany}, 10(2), 290-301.
#'
#'
#' #@export

prod_mod_ts <- function(data, method = "Schaefer", SSE_method = "standard",
                        B0_init = NA, B0_est = NA,
                        effort_unit = 1, plot = TRUE){

  res <- data
  year <- res$year
  Y <- res$Y
  if("f" %in% names(res)){
    f <- res$f
  }else f <- res$CPUE / Y
  if("CPUE"  %in% names(res)){
    CPUE <- res$CPUE
  }else CPUE <- Y/f

  # Initialisation of parameters (K, B0, r with fixed q)
  if(is.na(B0_init)){
    if(is.na(B0_est)){
      B0 <- 2 * mean(Y)
    }else B0 <- 2 * mean(Y[1:B0_est])
  }else B0 <- B0_init
  K <- B0 * 1.3
  r <- 1
  q <- mean(CPUE)/B0

  input <- c(K, B0, r)
  B <- B0

  # Evaluating model fit
  ssefn<-function(input){
    K <- input[1]
    B0 <- input[2]
    r <- input[3]

    B <- B0
    Yvec <- NULL
    Bvec <- NULL
    CPUE_hat <- NULL   # predicted CPUE
    yrs <- 1:length(Y)

    for (y in yrs){
      if(method == "Schaefer") SY <- r * B * (1 - B/K)
      if(method == "Fox") SY <- r * B * log(K/B)
      Bvec <- c(Bvec, B)
      CPUE_hat <- c(CPUE_hat, q * B)
      B <- B + SY - Y[y]
      B <- ifelse(B < 0, 0, B)
    }
    switch(SSE_method,
           "standard"={
             SSE <- sum((CPUE - CPUE_hat)^2, na.rm=TRUE)
           },
           "log"={
             SSE <- sum((log(CPUE) - log(CPUE_hat))^2, na.rm=TRUE)
           },
           "Thiels_U_statistic"={
             CPUE_diff <- rep(NA,length(yrs))
             for(i in 2:length(yrs)){
               CPUE_diff[i] <- CPUE[i] - CPUE[i-1]
             }
             SSE <- sqrt(sum((CPUE - CPUE_hat)^2, na.rm=TRUE) /
                               (sum((CPUE_diff)^2, na.rm=TRUE)))
           }
    )
    #SSE <- sum((CPUE - CPUE_hat)^2, na.rm=TRUE)  # take log? (with logarithm of parameters no unrealistic negative results are provided)
                                     # take Theil's U-statistic (Wittink 1988) ?? :
# ALternative SSE according to Theil's U-statistic:
#     CPUE_diff <- rep(NA,length(yrs))
#     for(i in 2:length(yrs)){
#       CPUE_diff[i] <- CPUE[i] - CPUE[i-1]
#     }
#     SSE <- sqrt( sum((CPUE - CPUE_hat)^2, na.rm=TRUE) /
#                   (sum((CPUE_diff)^2, na.rm=TRUE)))
    return(SSE)
  }

  estA <- nlm(ssefn, input, typsize = input, iterlim = 1000)
  estA <- nlm(ssefn, estA$est, typsize = input, iterlim = 1000)

  # Estimation of r and q
  K <- estA$est[1]
  B0 <- estA$est[2]
  r <- estA$est[3]
  q <- mean(CPUE) / B0

  input2 <- c(r, q)
  B <- B0

  ssefn2 <- function(input){
    r <- input[1]
    q <- input[2]
    B <- B0
    Yvec <- NULL
    Bvec <- NULL
    CPUE_hat <- NULL
    yrs <- 1:length(Y)
    for (y in yrs){
      if(method == "Schaefer") SY <- r * B * (1 - B/K)
      if(method == "Fox") SY <- r * B * log(K/B)
      Bvec <- c(Bvec, B)
      CPUE_hat <- c(CPUE_hat, q * B)
      B <- B + SY - Y[y]
      B <- ifelse(B<0,0,B)
    }
    switch(SSE_method,
           "standard"={
             SSE <- sum((CPUE - CPUE_hat)^2, na.rm=TRUE)
           },
           "log"={
             SSE <- sum((log(CPUE) - log(CPUE_hat))^2, na.rm=TRUE)
           },
           "Thiels_U_statistic"={
             CPUE_diff <- rep(NA,length(yrs))
             for(i in 2:length(yrs)){
               CPUE_diff[i] <- CPUE[i] - CPUE[i-1]
             }
             SSE <- sqrt(sum((CPUE - CPUE_hat)^2, na.rm=TRUE) /
                               (sum((CPUE_diff)^2, na.rm=TRUE)))
           }
    )
    #SSE <- sum((CPUE - CPUE_hat) ^ 2, na.rm = TRUE)
    # ALternative SSE according to Theil's U-statistic:
#     CPUE_diff <- rep(NA,length(yrs))
#     for(i in 2:length(yrs)){
#       CPUE_diff[i] <- CPUE[i] - CPUE[i-1]
#     }
#     SSE <- sqrt( sum((CPUE - CPUE_hat)^2, na.rm=TRUE) /
#                    (sum((CPUE_diff)^2, na.rm=TRUE)))
    return(SSE)
  }

  estB <- nlm(ssefn2, input2, typsize=input2, iterlim=1000)
  estB <- nlm(ssefn2, estB$est, typsize=estB$est, iterlim=1000)


  # calculate predicted biomass for all age classes
  K <- estA$est[1]
  B0 <- estA$est[2]
  r <- estB$est[1]
  q <- estB$est[2]
  B <- B0
  Yvec <- NULL
  Bvec <- NULL
  CPUE_hat <- NULL
  yrs <- 1:length(Y)
  for (y in yrs){
    if(method == "Schaefer") SY <- r * B * (1 - B/K)
    if(method == "Fox") SY <- r * B * log(K/B)
    Bvec <- c(Bvec, B)
    CPUE_hat <- c(CPUE_hat, q * B)
    B <- B + SY - Y[y]
  }

  # calculate MSY, Bmsy, Fmsy and Emsy
  if(method == "Schaefer"){
    MSY <- r * K / 4
    Bmsy <- K / 2
    Fmsy <- r / 2
    Emsy <- r / (2 * q) * effort_unit
  }
  if(method == "Fox"){
    MSY <- (K * r) / exp(1)
    Bmsy <- K / exp(1)
    Fmsy <- r / exp(1)
    Emsy <- r / (exp(1) * q) * effort_unit
  }


  # create ouput list
  ret <- list(
    year = year,
    Y = Y,
    f = f,
    method = method,
    CPUE = CPUE,
    Bvec = Bvec,
    CPUE_hat = CPUE_hat,
    K = K,
    r = r,
    q = q,
    MSY = MSY,
    Bmsy = Bmsy,
    Emsy = Emsy,
    Fmsy = Fmsy)
  class(ret) <- "prod_mod_ts"

  # create plot
  # if(plot) plot(ret)

  return(ret)
}


