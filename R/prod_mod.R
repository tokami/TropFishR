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
#' @param plot logical; if TRUE a graph is displayed
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
#'   \item \strong{Fox_Bv}: virgin biomass according to Fox model,
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
#'    to \link{prod_mod_dyn}.
#'
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
  }else f <- res$CPUE / Y
  if("CPUE"  %in% names(res)){
    CPUE <- res$CPUE
  }else CPUE <- Y/f

  mean.f <- mean(f,na.rm=T)
  sd.f <- sd(f,na.rm=T)

  # first: Schaefer, second: Fox
  CPUEs <- list(CPUE, log(CPUE))

  MSYs <- rep(NA,2)
  fMSYs <- rep(NA,2)
  as <- rep(NA,2)
  bs <- rep(NA,2)

  for(i in 1:2){
    CPUEi <- CPUEs[[i]]
    mean.CPUEi <- mean(CPUEi,na.rm=T)
    sd.CPUEi <- sd(CPUEi,na.rm=T)
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




