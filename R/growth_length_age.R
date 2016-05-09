#' @title Estimation from growth parameter using length-at-age data
#'
#' @description  This function estimates growth parameters from length-at-age data. It
#'    allows to perform different methods: Gulland and Holt, Ford Walford plot,
#'    Chapman's method, Bertalanffy plot, or non linear least squares method (LSM).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{age}: age measurements,
#'   \item \strong{length}: corresponding lengths.
#' }
#' @param method indicating which of following methods should be applied:
#'    \code{"GullandHolt"},
#'    \code{"FordWalford"}, \code{"Chapman"}, \code{"BertalanffyPlot"},
#'    or \code{"LSM"}
#' @param Linf_est BertalanffyPlot requires an estimate for Linf to derive K and t0
#'    (for more information see Details).
#' @param Linf_init initital parameter of Linf for non-linear sqaures fitting (default 10)
#' @param K_init initital parameter of K for non-linear sqaures fitting (default 0.1)
#' @param t0_init initital parameter of t0 for non-linear sqaures fitting (default 0)
#'
#' @examples
#' # synthetical length at age data
#' dat <- list(age = 1:7, length = c(25.7,36,42.9,47.5,50.7,52.8,54.2))
#' growth_length_age(dat, method = "GullandHolt")
#'
#' # Bertalaffy plot
#' dat <- list(age = c(0.64,1.16,1.65,2.1,2.64,3.21),
#'    length = c(17.3,27.9,35.3,40.2,43.3,45.5))
#' growth_length_age(dat, method = "BertalanffyPlot", Linf_est = 50)
#'
#' growth_length_age(dat, method = "LSM")
#'
#' @details
#' Gulland and Holt plot assumes
#' infinitestimal delta t (only reasonable approximation of growth parameters if delta t
#' is small). Ford Walford plot and Chapman assume constant time intervals between ages
#' (delta t). The Bertalanffy plot is a robust method, however it requires an estimate of Linf. As
#' long as this estimate is reasonable the resulting estimate of K is reasonable. For
#' a first estimate of Linf the Powell Wetherall method \link{powell_wetherall} can
#' be used. Otherwise, the largest fish or the average of the ten largest fish can be
#' used for a small or large sample, respectively. All lengths have to be smaller than
#' Linf as otherwise the logarithm is not defined. Oldest fish (if larger than Linf) have
#' to be omitted. Non-linear least squares fitting is the preferred method to estimate
#' growth parameters according to Sparre and Venema (1998).
#'
#' @return A list with the input parameters and following parameters:
#' \itemize{
#'  \item \strong{x}: independent variable used for regression analysis,
#'  \item \strong{y}: dependent variable used for regression analysis,
#'  \item \strong{reg_coeffs}: regression coefficients,
#'  \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'  \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'  \item \strong{t0}: theoretical time zero, at which individuals of this species hatch
#'      (only if Bertalanffy plot was applied).
#' }
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

growth_length_age <- function(param, method, Linf_est = NA,
                              Linf_init = 10, K_init = 0.1, t0_init = 0){

  res <- param
  t <- res$age
  Lt <- res$length

  delta_t <- diff(t)

  switch(method,
         "GullandHolt"={
           # growth rate
           y <- rep(NA,length(delta_t))
           for(i in 1:length(delta_t)){
             y[i] <- (Lt[i+1] - Lt[i])/delta_t[i]
           }

           #mean Lt
           x <- rep(NA,(length(Lt)-1))
           for(i in 1:(length(Lt)-1)){
             x[i] <- (Lt[i] + Lt[i+1])/ 2
           }

           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           K <- -b
           Linf <- -a/b

           plot(y ~ x, pch = 16,
                ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlab = "mean L(t)",
                ylab=expression(paste(delta,"L / ",delta,"t")),
                main = "Gulland and Holt")
           abline(a,b)

           tmp <- list(x=x,y=y)

         },
         "FordWalford"={
           if(!all(delta_t == 1)) stop("The Ford Walford method assumes constant time intervals!")

           delta_t <- delta_t[1]

           y <- Lt[2:length(Lt)]
           x <- Lt[-length(Lt)]

           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           K <- -1/delta_t * log(b)
           Linf <- a/(1-b)

           plot(y ~ x, pch = 16, ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(0,max(x,na.rm=TRUE)*1.1),
                xlab = "L(t)",
                ylab=expression(paste("L(t+",delta,"t)")),
                main = "Ford Walford plot")
           abline(a,b)
           abline(1,1, lty=2)
           segments(x0 = 0, y0 = Linf, x1 = Linf, y1 = Linf, lty=2)
           segments(x0 = Linf, y0 = 0, x1 = Linf, y1 = Linf, lty=2)

           tmp <- list(x=x,y=y)
         },
         "Chapman"={
           if(!all(delta_t == 1)) stop("The Chapman method assumes constant time intervals!")

           delta_t <- delta_t[1]

           Lt_short <- Lt[-length(Lt)]
           y <- Lt[2:length(Lt)] - Lt_short
           x <- Lt_short


           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           c <- -b
           K <- -(1/delta_t) * log(1+b)
           Linf <- -a/b

           plot(y ~ x, pch = 16, ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(0,max(x,na.rm=TRUE)*1.1),
                xlab = "L(t)",
                ylab=expression(paste("L(t+",delta,"t)-L(t)")),
                main = "Chapman")
           abline(a,b)

           tmp <- list(x=x,y=y)
         },
         "BertalanffyPlot"={

           #Linf should be given: otherwise using optim?
           if(is.na(Linf_est)) stop("For the Bertalanffy plot you have to provide an estimate for Linf_est!")
           y <- -log(1 - Lt / Linf_est)
           x <- t

           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           K <- b
           t0 <- a/-K
           Linf <- Linf_est

           plot(y ~ x, pch = 16, ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(0,max(x,na.rm=TRUE)*1.1),
                xlab = "t",
                ylab=expression(paste("-ln(1-L(t)/Linf)")),
                main = "Bertalanffy plot")
           abline(a,b)

           tmp <- list(x=x,y=y)


           if(max(Lt, na.rm=TRUE) > Linf_est) writeLines("Some lengths were larger as the Linf value which was \nprovided. These lengths are omitted.")
         },
         "LSM"={

           nls_mod <- nls(Lt ~ (Linf * (1 - exp(-K * (t - t0)))),
                      start = list(Linf = Linf_init, K = K_init, t0 = t0_init))

           plot(Lt ~ t,
                ylab = "L(t)",
                xlab = "t(age)",
                main = "Non-linear least squares method")
           lines(t, predict(nls_mod))
           sum_mod <- summary(nls_mod)
           Linf <- sum_mod$coefficients[1]
           K <- sum_mod$coefficients[2]
           t0 <- sum_mod$coefficients[3]

           tmp <- list(t0 =t0)
         })

  res2 <- c(res,list(delta_t = delta_t))
  if(exists("x", where = tmp)) res2$x <- x
  if(exists("y", where = tmp)) res2$y <- y

  ret <- c(res2,list(reg_coeffs = sum_mod$coefficients,
                    Linf = Linf,
                    K = K))

  if(exists("t0", where = tmp)) ret$t0 <- t0
  return(ret)
}

