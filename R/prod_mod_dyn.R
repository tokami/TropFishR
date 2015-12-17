#' @title Non-equilibrium dynamic production models
#'
#' @description Non-equilibrium dynamic biomass models
#'
#' @param data a dataframe of parameters
#' \itemize{
#'   \item \strong{year} years,
#'   \item \strong{yield} catch in weight of fishery per year, and
#'   \item \strong{effort} fishing effort per year, or
#'   \item \strong{CPUE} catch per unit of effort per year.
#' }
#'
#' @keywords function biomass MSY stock
#'
#' @examples
#' # load data
#' data(emperor)
#'
#' # run model
#' prod_mod_dyn(data = emperor)
#'
#' @details Either the Catch per unit of effort (CPUE) is inserted
#'    into the model directly (by a column \strong{CPUE}) or the CPUE is calculated from
#'    the catch and effort, then these two vectors should have required units.
#'
#' @references
#' Dharmendra, D., Sólmundsson, J., 2005. Stock assessment of the offshore Mauritian banks
#' using dynamic biomass models and analysis of length frequency of the Sky Emperor
#' (\emph{Lethrinus mahsena}). Fisheries Training Program The United Nations University, 61
#'
#' Hilborn, R. and Walters, C., 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York
#'
#'
#' @export

prod_mod_dyn <- function(data){

  res <- data
  year <- res$year
  Y <- res$Y
  if("f" %in% names(res)){
    f <- res$f
  }else f <- res$CPUE / Y
  if("CPUE"  %in% names(res)){
    CPUE <- res$CPUE
  }else CPUE <- Y/f

  # FORMULA== By+1 = By + r By(1 - By/K) - Yy
  # Yy=yield/catch
  # Initial parameters of starting biomass B0, Carrying capacity K and the rate of pop. growth r
  # q= coefficient catchability
  # There are four parameters to be estimated K, B0,r and q
  #It is not advisable to estimate all- first estimate K, B0 and r with fixed q, then estimate r and q, then all
  B0 <- 2*mean(Y)
  K <- B0*1.3
  r <- 1
  q <- mean(CPUE)/B0

  # biomass greater than catch # K>B0- carrying capacity
  # rate of population growth # as I=qB
  input <- c(K,B0,r)
  B <- B0

  ###
  ssefn<-function(input){
    K <- input[1]
    B0 <- input[2]
    r <- input[3]

    B <- B0
    Yvec <- NULL
    Bvec <- NULL
    Ihat <- NULL
    yrs <- 1:length(Y)

    for (y in yrs){
      SY <- r * B * (1 - B/K)
      Bvec <- c(Bvec, B)
      Ihat <- c(Ihat, q * B)
      B <- B + SY - Y[y]
      B <- ifelse(B < 0, 0, B)
    }
    SSE <- sum((CPUE - Ihat)^2)
    return(SSE)
  }

  # to optimise
  estA <- nlm(ssefn, input, typsize=input, iterlim=1000)
  estA <- nlm(ssefn, estA$est, typsize=input, iterlim=1000)
  # nlm- non linear minimization
  # using the result of the first estimate we estimate again
  # estimate r and q
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
    # Estimation of the two parameters:
    Bvec <- NULL
    Ihat <- NULL
    yrs <- 1:length(Y)
    for (y in yrs){
      SY <- r * B * (1 - B/K)
      Bvec <- c(Bvec, B)
      Ihat <- c(Ihat, q * B)
      B <- B + SY - Y[y]
      B <- ifelse(B<0,0,B)
    }
    SSE <- sum((CPUE - Ihat) ^ 2)
    # really necessary?
     #plot(yrs,CPUE, type="b", xlab="year", ylab="CPUE")
     #lines(Ihat, col="red")
    return(SSE)
  }

  estB <- nlm(ssefn2, input2, typsize=input2, iterlim=1000)
  estB <- nlm(ssefn2, estB$est, typsize=estB$est, iterlim=1000)

  #####
  # plotting
  # to calculate the predicted biomass for all ages and the predicted index
  K <- estA$est[1]
  B0 <- estA$est[2]
  r <- estB$est[1]
  q <- estB$est[2]
  B <- B0
  Yvec <- NULL
  Bvec <- NULL
  Ihat <- NULL
  yrs <- 1:length(Y)
  for (y in yrs){
    SY <- r * B * (1 - B/K)
    Bvec <- c(Bvec, B)
    Ihat <- c(Ihat, q * B)
    B <- B + SY - Y[y]
  }

  # calculate MSY, Bmsy, Fmsy and Emsy
  MSY <- r*K/4
  Bmsy <- K/2
  Emsy <- r/(2*q)*1000  #######   why 1000 ?????
  Fmsy <- r/2

  # create ouput list
  ret <- list(
    year = year,
    Y = Y,
    f = f,
    CPUE = CPUE,
    Bvec = Bvec,
    Ihat = Ihat,
    K = K,
    r = r,
    q = q,
    MSY = MSY,
    Bmsy = Bmsy,
    Emsy = Emsy,
    Fmsy = Fmsy)
  class(ret) <- "prod_mod_dyn"

  # create plot
  plot(ret)

  return(ret)
}


