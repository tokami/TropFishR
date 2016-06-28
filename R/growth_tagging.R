#' @title Growth from tagging data
#'
#' @description  This function estimates growth parameters from tagging data. Munro plot
#'    is applied
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{L1}: length at tagging [cm],
#'   \item \strong{L2}: length at recapture [cm],
#'   \item \strong{delta_t}: time interval between tagging ang recapture
#'   (instead two vectors with \strong{t1} (age at tagging) and \strong{t2}
#'   (age at recapture) can be provided.
#' }
#' @param method indicating which of following methods should be applied: "GullandHolt"
#'    or "Munro".
#' @param Linf_range two values indicating the lower and upper limits of the range,
#'    in which the \link{optimise} searches for the Linf value with the best fit
#'    (lowest CV value ),
#' @param time_unit indicating the unit of the time interval, either "year", "month",
#'    "week", or "day"
#'
#' @examples
#' # from Wolff (1984)
#' dat <- list(L1 = c(40,46,29,30,18,31,48,49,59,58,61,65,57,55),
#'    L2 = c(85,53,55,56,25,43,70,59,62,80,72,83,65,56),
#'    delta_t = c(289,26,84,77,14,38,89,38,28,149,89,74,38,21))
#' growth_tagging(param = dat, "Munro", time_unit = "day", Linf_range=c(80,120))
#' growth_tagging(param = dat, "GullandHolt", time_unit = "day")
#'
#'
#' # from Sparre and Venema (1999)
#' dat <- list(L1 = c(9.7,10.5,10.9,11.1,12.4,12.8,14.0,16.1,16.3,17.0,17.7),
#'    L2 = c(10.2,10.9,11.8,12.0,15.5,13.6,14.3,16.4,16.5,17.2,18.0),
#'    delta_t = c(53,33,108,102,272,48,53,73,63,106,111))
#' growth_tagging(param = dat, "Munro", time_unit = "day", Linf_range = c(10,40))
#' growth_tagging(param = dat, "GullandHolt", time_unit = "day")
#'
#' @details
#' If Munro plot is applied the optimal Linf value is found by minimizing the coefficient
#' of variation (CV = sd(K)/mean(K)). For this iterative method the \link{optimise}
#' function is applied. The histogram of the individual K values allows to distinguish
#' potential differences in growth performance between individuals. t0 can not be
#' estimated by Munro plot, neither by the Gulland Holt method.
#'
#' @return A list with the input parameters and following parameters:
#' \itemize{
#'  \item \strong{x}: independent variable used for regression analysis,
#'  \item \strong{y}: dependent variable used for regression analysis,
#'  \item \strong{reg_coeffs}: regression coefficients,
#'  \item \strong{r2}: r squared of regression analysis,
#'  \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'  \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'  \item \strong{conf_int_K}: confidence intervals of K (only if Gulland Holt method
#'      was applied).
#' }
#'
#' @importFrom grDevices colorRampPalette dev.new dev.off recordPlot rgb
#' @importFrom graphics abline axis barplot box contour grid hist identify image layout legend lines locator matplot mtext par plot points rect segments text title
#' @importFrom stats aggregate deviance dnorm lm lowess na.omit nlm nls optim optimise optimize predict qt rnorm sd update
#'
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock assessment.
#' Part 2. Excercises. FAO Fisheries Technical Paper, (306.2, Rev. 2). 94 p.
#'
#' Wolff, M., 1984. Early setback for scallop culture in Peru.
#'
#' @export

growth_tagging <- function(param, method,
                            Linf_range = c(5,600),time_unit = "year"){

  res <- param

  L1 <- res$L1
  L2 <- res$L2
  if("delta_t" %in% names(res)) delta_t <- res$delta_t
  if(!"delta_t" %in% names(res)) delta_t <- res$t2 - res$t1

  switch(time_unit,
         "year" ={
           time_u <- 1
         },
         "month"={
           time_u <- 12
         },
         "week"={
           time_u <- 52
         },
         "day"={
           time_u <- 365
         })

  # growth increment per time / growth rate
  incr_time <- ((L2-L1) / delta_t ) * time_u

  switch(method,
         "GullandHolt"={
           y <- incr_time

           #mean Lt
           x <- (L1 + L2) / 2

           sum_mod <- summary(lm(y ~ x))
           coeffs <- sum_mod$coefficients
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           K <- -b
           Linf <- a/K
           R2 <- sum_mod$r.squared
           # standard error and confidence limits
           sb2 <- (1/(length(L1)-2)) * ((sd(y,na.rm = TRUE)/sd(x,na.rm = TRUE))^2 - b^2)
           sb <- sqrt(sb2)
           SE_b <- abs(b) * qt(0.975,sum_mod$df[2])
           tg <- qt(0.975,sum_mod$df[2])
           conf_int <- c(K - tg * sb,K + tg * sb)

           plot(y ~ x, pch = 16,
                ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(min(x,na.rm=TRUE)*0.9,Linf*1.05),
                xlab = "mean L(t)",
                ylab=expression(paste(delta,"L / ",delta,"t")),
                main = "Gulland and Holt")
           abline(a,b)

           tmp <- list(y = y, x = x, coeffs = coeffs, R2 = R2, conf_int = conf_int)
         },

         "Munro"={
           func <- function(Linf){
             K <- ((log(Linf-L1) - log(Linf-L2)) / (delta_t)) * time_u
             CV <- sd(K, na.rm = TRUE) / mean(K, na.rm = TRUE)
           }

           opt_mod <- optimise(f = func, interval = Linf_range)
           Linf <- opt_mod$minimum
           K <- opt_mod$objective
           indi_Ks <- ((log(Linf-L1) - log(Linf-L2)) / (delta_t)) * time_u
           hist(indi_Ks, breaks = length(L1),
                xlab = "K",
                main = "Histogram  of individual Ks")

           tmp <- list(indi_Ks = indi_Ks)
         })

  res2 <- res
  if(exists("x", where = tmp)) res2$x <- x
  if(exists("y", where = tmp)) res2$y <- y
  if(exists("coeffs", where = tmp)) res2$reg_coeffs <- coeffs
  if(exists("R2", where = tmp)) res2$r2 <- R2
  if(exists("indi_Ks", where = tmp)) res2$indi_Ks <- indi_Ks
  ret <- c(res2,list(Linf = Linf, K = K))
  if(exists("conf_int", where = tmp)) ret$conf_int_K <- conf_int
  return(ret)
}
