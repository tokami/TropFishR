#' @title Production models
#'
#' @description Production models are holistic models, which
#'    can be used to estimate maximum sustainable yield (MSY) and virgin biomass. This
#'    function uses the equilibrium approach to estimate parameters (Schaefer model and
#'    Fox model).
#'
#' @param data a dataframe consisting of:
#' \itemize{
#'   \item \code{year} year vector,
#'   \item \code{Y} catch in weight per year, and
#'   \item \code{f} fishing effort per year, or
#'   \item \code{CPUE} catch per unit of effort per year (optional).
#' }
#' @param plot logical; if TRUE, a graph is displayed
#'
#' @keywords function biomass MSY production surplus
#'
#' @examples
#' data(trawl_fishery_Java)
#' prod_mod(data = trawl_fishery_Java, plot = TRUE)
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{Schaefer_lm}: intercept and slope of linear model following
#'   the Schaefer model,
#'   \item \strong{Fox_lm}: intercept and slope of linear model following
#'   the Fox model,
#'   \item \strong{Schaefer_MSY}: MSY according to Schaefer model,
#'   \item \strong{Schaefer_fMSY}: fishing effort yielding in MSY according to
#'   Schaefer model,
#'   \item \strong{Schaefer_Bv}: virgin biomass according to Schaefer model,
#'   \item \strong{ln_CPUE}: natural logarithm of CPUE values,
#'   \item \strong{Fox_MSY}: MSY according to Fox model,
#'   \item \strong{Fox_fMSY}: fishing effort yielding in MSY according to Fox model,
#'   \item \strong{Fox_Bv}: virgin biomass according to Fox model.
#' }
#'
#' @details Production models are also called surplus production models or
#'    biomass dynamic models. They can be applied if sufficient data are available: effort
#'    and yield parameters have to be expended over a certain number of years. Furthermore,
#'    the fishing effort must have undergone substantial changes over the period covered
#'    (Sparre and Venema, 1998). Either the catch per unit of effort (CPUE) is inserted
#'    into the model directly (objectname: \code{CPUE}) or the CPUE is calculated from
#'    the catch and effort, then these two vectors should have required units. There
#'    are three ways of estimating paramaters of production models, (i) assuming
#'    equlibrium conditions, (ii) transforming equation to linear form, or (iii)
#'    time-series fitting (Hilborn and Walters, 1992). The first approach corresponds
#'    to the Schaefer and Fox model and thus the methodology of this function.
#'    The authors recommend to use dynamic fitting methods when possible rather than
#'    the equilibrium approach. For dynamic production models please refer
#'    to \link{prod_mod_ts}.
#'
#' @importFrom graphics plot
#' @importFrom stats lm qt sd
#'
#' @references
#' Fox, W. W. Jr., 1970. An exponential surplus-yield model for optimizing exploited fish
#' populations. \emph{Trans.Am.Fish.Soc.}, 99:80-88
#'
#' Graham, M., 1935. Modern theory of exploiting a fishery and application to North Sea
#' trawling. \emph{J.Cons.CIEM}, 10(3):264-274
#'
#' Hilborn, R., Walters, C. J. (1992). Quantitative fisheries stock assessment:
#' choice, dynamics and uncertainty. \emph{Reviews in Fish Biology and Fisheries}, 2(2),
#' 177-178.
#'
#' Schaefer, M., 1954. Some aspects of the dynamics of populations important to the
#' management of the commercial marine fisheries. \emph{Bull.I-ATTC/Bol. CIAT}, 1(2):27-56
#'
#' Schaefer, M., 1957. A study of the dynamics of the fishery for yellowfin tuna of the
#' eastern tropical Pacific Ocean [in English and Spanish]. \emph{Ibid.}, 2(6): 245-285
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

prod_mod <- function(data, plot = FALSE){
  res <- data
  year <- res$year
  Y <- res$Y
  if("f" %in% names(res)){
    f <- res$f
  }else f <- Y/res$CPUE
  if("CPUE"  %in% names(res)){
    CPUE <- res$CPUE
  }else CPUE <- Y/f

  mean.f <- mean(f,na.rm=TRUE)
  sd.f <- sd(f,na.rm=TRUE)

  # first: Schaefer, second: Fox
  CPUEs <- list(CPUE, log(CPUE))

  MSYs <- rep(NA,2)
  fMSYs <- rep(NA,2)
  as <- rep(NA,2)
  bs <- rep(NA,2)

  for(i in 1:2){
    CPUEi <- CPUEs[[i]]
    mean.CPUEi <- mean(CPUEi,na.rm=TRUE)
    sd.CPUEi <- sd(CPUEi,na.rm=TRUE)
    mod.sum <- summary(lm(CPUEi ~ f))
    a <- mod.sum$coefficients[1]
    b <- mod.sum$coefficients[2]
    sb2 <- ((sd.CPUEi/sd.f)^2 - b^2 )/ (length(f)-2)
    sb <- mod.sum$coefficients[4]
    tdis <- qt(0.975, mod.sum$df[2])
    conf.b <- c(b - tdis * sb, b + tdis * sb)
    sa2 <- sb2* (sd.f^2 * (length(f)-1)/length(f) + mean.f^2)
    sa <- mod.sum$coefficients[3]
    conf.a <- c(a- tdis * sa, a+ tdis * sa)

    if(i == 1){
      MSY <- -0.25 * a^2/b
      fMSY <- -0.5 * a/b
    }
    if(i == 2){
      MSY <- -(1/b) * exp(a-1)
      fMSY <- -1/b
    }
    as[i] <- a
    bs[i] <- b
    MSYs[i] <- MSY
    fMSYs[i] <- fMSY
  }

  ret <- list(
    year = year,
    Y = Y,
    f = f,
    CPUE = CPUEs[[1]],
    Schaefer_lm = c(as[1],bs[1]),
    Fox_lm = c(as[2],bs[2]),
    Schaefer_MSY = MSYs[1],
    Schaefer_fMSY = fMSYs[1],
    Schaefer_Bv = as[1],
    ln_CPUE = CPUEs[[2]],
    Fox_MSY = MSYs[2],
    Fox_fMSY = fMSYs[2],
    Fox_Bv = as[2])

  class(ret) <- "prod_mod"

  #Plot
  if(plot) plot(ret)

  return(ret)
}





#' @title  Production models with time series fitting
#'
#' @description This function applies the production models under non-equilibrium
#'      conditions by applying time series fitting using non-linear least squares
#'      minimisation.
#'
#' @param data a dataframe of parameters
#' \itemize{
#'   \item \code{year} years,
#'   \item \code{yield} catch in weight of fishery per year,
#'   \item \code{effort} fishing effort per year,
#'   \item \code{CPUE} catch per unit of effort per year (optional).
#' }
#' @param method indicating if Schaefer or Fox model should be applied. First assumes a
#'    logistic relationship between growth rate and biomass, whereas second assumes it to
#'    foolow the Gompertz distribution (Richards 1959). Default is the dynamic Schaefer model.
#' @param B0_init numeric; if realistic initial estimate for virgin biomass is available.
#'    If NA initial estimate for virgin biomass is set to two times average yield of all
#'    or part of yield values (see \code{B0_est}).
#' @param B0_est intital value of virgin biomass estimating using all yield values (NA) or
#'    first years of time series, then provide numerical representing number of years
#' @param effort_unit multiplication factor for the unit of effort. Default is 1.
#' @param plot logical; if TRUE (default) a graph is displayed
#'
#' @keywords function biomass MSY production
#'
#' @examples
#' data(emperor)
#' prod_mod_ts(emperor, method = "Schaefer")
#' prod_mod_ts(emperor, method = "Fox")
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{Bvec}: biomass vector,
#'   \item \strong{CPUE_hat}: predicted CPUE,
#'   \item \strong{K}: carrying capacity,
#'   \item \strong{r}: population growth rate,
#'   \item \strong{q}: catchability coefficient,
#'   \item \strong{MSY}: maximum sustainabale yield (MSY),
#'   \item \strong{Bmsy}: biomass at MSY,
#'   \item \strong{Emsy}: fishing effort at MSY
#'   \item \strong{Fmsy}: fishing mortality at MSY,
#' }
#'
#' @details Either catch per unit of effort (CPUE) is inserted
#'    into the model directly (by a column \code{CPUE}) or CPUE is calculated from
#'    the catch and effort, then these two vectors should have required units.
#'    Whenever a good estimate for the virigin biomass is available, this estimate
#'    should be inserted for \code{B_init}. The default approach for the initial
#'    estimate of the virgin biomass is to multiply the average yield by 2 (Dharmendra
#'    and Solmundsson, 2005). Alternatively, just a part of the time series of
#'    yield values can be choosen to represent the virgin biomass.
#'    The minimisation procedure is based on least error sum of squares (SSE). For
#'    the logistic (Schaefer) method the standard calculation of SSE is applied
#'    (\code{sum((CPUE - predicted CPUE)^2)}), for
#'    the method with Gompertz distribution (Fox) SSE is calculated according to
#'    the Thiel's U statistic \code{sqrt(sum(CPUE - predicted CPUE)/sum(CPUE(t) - CPUE(t-1)))}
#'    (Wittink, 1988).
#'
#' @importFrom graphics plot
#' @importFrom stats nlm
#'
#' @references
#' Dharmendra, D., Solmundsson, J., 2005. Stock assessment of the offshore Mauritian banks
#' using dynamic biomass models and analysis of length frequency of the Sky Emperor
#' (\emph{Lethrinus mahsena}). Fisheries Training Program The United Nations University, 61
#'
#' Hilborn, R. and Walters, C., 1992. Quantitative Fisheries Stock Assessment: Choice,
#' Dynamics and Uncertainty. Chapman and Hall, New York
#'
#' Prager, M. H., 1994. A suite of extensions to a non-equilibrium surplus production
#' model. \emph{Fishery Bulletin} 92: 374-389
#'
#' Richards, F. J., 1959. A flexible growth function for empirical use.
#' \emph{Journal of experimental Botany}, 10(2), 290-301.
#'
#' Wittink, D. R., 1988. The application of regression analysis. Allyn and Bacon. Inc.
#' Boston. MA. 324p.
#'
#'
#' @export

prod_mod_ts <- function(data, method = "Schaefer",
                        B0_init = NA, B0_est = NA,
                        effort_unit = 1, plot = TRUE){

  res <- data
  year <- res$year
  Y <- res$Y
  if("f" %in% names(res)){
    f <- res$f
  }else f <- Y/res$CPUE
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
    switch(method,
           "Schaefer"={
             SSE <- sum((CPUE - CPUE_hat)^2, na.rm=TRUE)
           },
           "Fox"={
             CPUE_diff <- rep(NA,length(yrs))
             for(i in 2:length(yrs)){
               CPUE_diff[i] <- CPUE[i] - CPUE[i-1]
             }
             SSE <- sqrt(sum((CPUE - CPUE_hat)^2, na.rm=TRUE) /
                               (sum((CPUE_diff)^2, na.rm=TRUE)))
           }
    )
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
    switch(method,
           "Schaefer"={
             SSE <- sum((CPUE - CPUE_hat)^2, na.rm=TRUE)
           },
           "Fox"={
             CPUE_diff <- rep(NA,length(yrs))
             for(i in 2:length(yrs)){
               CPUE_diff[i] <- CPUE[i] - CPUE[i-1]
             }
             SSE <- sqrt(sum((CPUE - CPUE_hat)^2, na.rm=TRUE) /
                               (sum((CPUE_diff)^2, na.rm=TRUE)))
           }
    )
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
  if(plot) plot(ret)

  return(ret)
}



