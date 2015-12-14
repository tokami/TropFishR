#' @title Production models - Schaefer and Fox models
#'
#' @description Surplus production models (or production models) are holistic models, which
#'    can be used to estimate the maximum sustainable yield, biomass and effort for a
#'    certain population. This function applies two models simultaneously: the Schaefer
#'    and the Fox model.
#'
#' @param data a dataframe with at least two columns:
#' \itemize{
#'   \item \strong{yield} catch in weight of fishery per year,
#'   \item \strong{effort} fishing effort per year,
#' }
#'
#' @examples
#' # load data
#' data(trawl_fishery_Java)
#'
#' # run model
#' prod_mod(data = trawl_fishery_Java)
#'
#' @details Production models can be applied if sufficient data are available: effort
#'    and yield parameters have to be expended over a certain number of years. Furthermore,
#'    the fishing effort must have undergone substantial changes over the period covered
#'    (Sparre and Venema, 1998).
#'
#' @references
#' Fox, W. W. Jr., 1970. An exponential surplus-yield model for optimizing exploited fish
#' populations. \emph{Trans.Am.Fish.Soc.}, 99:80-88
#'
#' Graham, M., 1935. Modern theory of exploiting a fishery and application to North Sea
#' trawling. \emph{J.Cons.CIEM}, 10(3):264-274
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

prod_mod <- function(data){
  res <- data
  Y <- res$yield
  f <- res$effort

  mean.f <- mean(f,na.rm=T)
  sd.f <- sd(f,na.rm=T)

  #Schaefer
  CPUE.S <- Y/f
  mean.CPUE.S <- mean(CPUE.S,na.rm=T)
  sd.CPUE.S <- sd(CPUE.S,na.rm=T)
  a.S <- summary(lm(CPUE.S ~ f))$coefficients[1]
  b.S <- summary(lm(CPUE.S ~ f))$coefficients[2]
  sb2.S <- ((sd.CPUE.S/sd.f)^2 - b.S^2 )/ (length(f)-2)
  sb.S <- summary(lm(CPUE.S ~ f))$coefficients[4]
  tdis <- qt(0.975,summary(lm(CPUE.S ~ f))$df[2])
  conf.b.S <- c(b.S - tdis * sb.S, b.S + tdis * sb.S)
  sa2.S <- sb2.S * (sd.f^2 * (length(f)-1)/length(f) + mean.f^2)
  sa.S <- summary(lm(CPUE.S ~ f))$coefficients[3]
  conf.a.S <- c(a.S - tdis * sa.S, a.S + tdis * sa.S)
  MSY.S <- -0.25 * a.S^2/b.S
  fMSY.S <- -0.5 * a.S/b.S

  #Fox
  CPUE.F <- log(Y/f)
  mean.CPUE.F <- mean(CPUE.F,na.rm=T)
  sd.CPUE.F <- sd(CPUE.F,na.rm=T)
  a.F <- summary(lm(CPUE.F ~ f))$coefficients[1]
  b.F <- summary(lm(CPUE.F ~ f))$coefficients[2]
  sb2.F <- ((sd.CPUE.F/sd.f)^2 - b.F^2 )/ (length(f)-2)
  sb.F <- summary(lm(CPUE.F ~ f))$coefficients[4]
  conf.b.F <- c(b.F - tdis * sb.F, b.F + tdis * sb.F)
  sa2.F <- sb2.F * (sd.f^2 * (length(f)-1)/length(f) + mean.f^2)
  sa.F <- summary(lm(CPUE.F ~ f))$coefficients[3]
  conf.a.F <- c(a.F - tdis * sa.F, a.F + tdis * sa.F)
  MSY.F <- -(1/b.F) * exp(a.F-1)
  fMSY.F <- -1/b.F

  ret <- c(res,list(
    Schaefer_lm = c(a.S,b.S),
    Fox_lm = c(a.F,b.F),
    Schaefer_CPUE = CPUE.S,
    Schaefer_MSY = MSY.S,
    Schaefer_fMSY = fMSY.S,
    Fox_MSY = MSY.F,
    Fox_fMSY = fMSY.F
  ))
  class(ret) <- "prod_mod"

  #Plot
  plot(ret)

  return(ret)
}




