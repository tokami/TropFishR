#' @title Estimation of growth parameter using length-at-age data
#'
#' @description  This function estimates growth parameters from
#'    length-at-age data. It
#'    allows to perform different methods: Gulland and Holt, Ford Walford plot,
#'    Chapman's method, Bertalanffy plot, or non linear least squares method (LSM).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{age}: age measurements,
#'   \item \strong{length}: corresponding lengths in cm.
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
#' @param CI logical; Should confidence intervals be calculated? This option only
#'    works for the LSM method. Default is FALSE.
#' @param ci.level required confidence level (for LSM method only)
#' @param age_plot sequence with ages used for plotting (LSM method only). By default
#'      age_plot = seq(min(param$age),max(param$age),0.1)
#' @param do.sim logical. Should Monte Carlo simulation be applied? Default = FALSE
#' @param nsim the number of Monte Carlo simulations to be performed,
#'    minimum is 10000 (default).
#'
#' @examples
#' # synthetical length at age data
#' dat <- list(age = rep(1:7,each = 5),
#'    length = c(rnorm(5,25.7,0.9),rnorm(5,36,1.2),rnorm(5,42.9,1.5),rnorm(5,47.5,2),
#'    rnorm(5,50.7,0.4),rnorm(5,52.8,0.5),rnorm(5,54.2,0.7)))
#' growth_length_age(dat, method = "GullandHolt")
#'
#' # Bertalaffy plot
#' growth_length_age(dat, method = "BertalanffyPlot", Linf_est = 50)
#'
#' # non linear least squares method
#' \donttest{
#' output <- growth_length_age(param = dat, method = "LSM",
#'      Linf_init = 30, CI = TRUE, age_plot=NULL)
#' summary(output$mod)
#' }
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
#' growth parameters according to Sparre and Venema (1998). If \code{CI = TRUE} the
#' confidence interval of parameters is calculated and plotted. For plotting the
#' confidence interval the \code{\link{predictNLS}} from the \link{propagate} package
#' is applied.
#'
#' @return A list with the input parameters and following parameters:
#' \itemize{
#'  \item \strong{x}: independent variable used for regression analysis,
#'  \item \strong{y}: dependent variable used for regression analysis,
#'  \item \strong{mod}: (non) linear model,
#'  \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'  \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'  \item \strong{t0}: theoretical time zero, at which individuals of this species hatch
#'      (only for Bertalanffy plot and LSM method).
#'  \item \strong{estimates}: dataframe with growth parameters and confidence intervals
#'      (only if LSM method was applied).
#' }
#'
#' @importFrom propagate predictNLS
#' @importFrom graphics abline lines plot segments
#' @importFrom stats lm nls predict aggregate confint
#' @importFrom grDevices adjustcolor
#' @import MASS
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

growth_length_age <- function(param, method, Linf_est = NA,
                              Linf_init = 10, K_init = 0.1, t0_init = 0,
                              CI = FALSE, ci.level = 0.95,
                              age_plot = NULL,
                              do.sim = FALSE,
                              nsim = 10000){

  res <- param
  t <- res$age
  Lt <- res$length

  switch(method,
         "GullandHolt"={
           # Can not deal with multiple observations (lengths) per age
           if(length(unique(t)) == length(t)){
             delta_t <- diff(t)
           }else{
             delta_t <- diff(unique(t))
             Lt <- aggregate(Lt,list(t = t),mean,na.rm= TRUE)$x
           }

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

           mod <- lm(y ~ x)
           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           R2 <- sum_mod$r.squared
           K <- -b
           Linf <- -a/b

           plot(y ~ x, pch = 16,
                ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlab = "mean L(t)",
                ylab=expression(paste(delta,"L / ",delta,"t")),
                main = "Gulland and Holt")
           abline(a,b)

           tmp <- list(delta_t=delta_t,x=x,y=y,R2=R2)

         },
         "FordWalford"={
           # Can not deal with multiple observations (lengths) per age
           if(length(unique(t)) == length(t)){
             delta_t <- diff(t)
           }else{
             delta_t <- diff(unique(t))
             Lt <- aggregate(Lt,list(t = t),mean,na.rm= TRUE)$x
           }

           if(!all(delta_t == 1)) stop("The Ford Walford method assumes constant time intervals!")

           delta_t <- delta_t[1]

           y <- Lt[2:length(Lt)]
           x <- Lt[-length(Lt)]

           mod <- lm(y ~ x)
           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           R2 <- sum_mod$r.squared
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

           tmp <- list(delta_t=delta_t,x=x,y=y,R2=R2)
         },
         "Chapman"={
           # Can not deal with multiple observations (lengths) per age
           if(length(unique(t)) == length(t)){
             delta_t <- diff(t)
           }else{
             delta_t <- diff(unique(t))
             Lt <- aggregate(Lt,list(t = t),mean,na.rm= TRUE)$x
           }

           if(!all(delta_t == 1)) stop("The Chapman method assumes constant time intervals!")

           delta_t <- delta_t[1]

           Lt_short <- Lt[-length(Lt)]
           y <- Lt[2:length(Lt)] - Lt_short
           x <- Lt_short


           mod <- lm(y ~ x)
           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           R2 <- sum_mod$r.squared
           c <- -b
           K <- -(1/delta_t) * log(1+b)
           Linf <- -a/b

           plot(y ~ x, pch = 16, ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(0,max(x,na.rm=TRUE)*1.1),
                xlab = "L(t)",
                ylab=expression(paste("L(t+",delta,"t)-L(t)")),
                main = "Chapman")
           abline(a,b)

           tmp <- list(delta_t=delta_t,x=x,y=y,R2=R2)
         },
         "BertalanffyPlot"={

           delta_t <- diff(t)

           #Linf should be given: otherwise using optim?
           if(is.na(Linf_est)) stop("For the Bertalanffy plot you have to provide an estimate for Linf_est!")
           y <- -log(1 - Lt / Linf_est)
           x <- t

           mod <- lm(y ~ x)
           sum_mod <- summary(lm(y ~ x))
           a <- sum_mod$coefficients[1]
           b <- sum_mod$coefficients[2]
           R2 <- sum_mod$r.squared
           K <- b
           t0 <- a/-K
           Linf <- Linf_est

           plot(y ~ x, pch = 16, ylim = c(0,max(y,na.rm=TRUE)*1.1),
                xlim = c(0,max(x,na.rm=TRUE)*1.1),
                xlab = "t",
                ylab=expression(paste("-ln(1-L(t)/Linf)")),
                main = "Bertalanffy plot")
           abline(a,b)

           tmp <- list(x=x,y=y,R2=R2,t0=t0)


           if(max(Lt, na.rm=TRUE) > Linf_est) writeLines("Some lengths were larger as the Linf value which was \nprovided. These lengths are omitted.")
         },
         "LSM"={

           nls_mod <- nls(Lt ~ (Linf * (1 - exp(-K * (t - t0)))),
                      start = list(Linf = Linf_init, K = K_init, t0 = t0_init))

           if(!CI){
             plot(Lt ~ t,
                  ylab = "L(t)",
                  xlab = "t(age)",
                  main = "Non-linear least squares method")
             lines(t, predict(nls_mod))
           }
           sum_mod <- summary(nls_mod)
           Linf <- sum_mod$coefficients[1]
           K <- sum_mod$coefficients[2]
           t0 <- sum_mod$coefficients[3]

           mod <- nls_mod
           tmp <- list(t0 = t0)

           if(CI){
             suppressMessages(cis <- confint(nls_mod, level = ci.level))
             nls_res <- data.frame(Names=c("Linf","K","t0"),
                                   Value=c(round(Linf,2),
                                           round(K,2),
                                           round(t0,2)),
                                   Lower.CI = c(round(cis[1,1],2),
                                                round(cis[2,1],2),
                                                round(cis[3,1],2)),
                                   Upper.CI = c(round(cis[1,2],2),
                                                round(cis[2,2],2),
                                                round(cis[3,2],2)))
             tmp$nls_res = nls_res

             # plot with confidence interval
             if(is.null(age_plot)){
               age_plot <- seq(min(floor(t)),max(ceiling(t)),0.1)
             }else age_plot <- age_plot
             ## Taylor error propagation and Monte Carlo simulation for confidence interval
             ## if (!requireNamespace("propagate", quietly = TRUE)) {
             ##     stop("Package \"propagate\" needed for this function to work. Please install it.",
             ##          call. = FALSE)
             ## }
             sink(tempfile())
             pred_L <- suppressMessages(propagate::predictNLS(nls_mod,
                                                              do.sim = do.sim,
                                                              nsim = nsim,
                                             newdata = data.frame(t = age_plot)))
             sink()
             # Taylor propagation
             if(!do.sim){
               predVals <- cbind(age_plot,as.data.frame(pred_L$summary[,c(1,5,6)]))
             }
             # Monte Carlo simulation
             if(do.sim){
               predVals <- cbind(age_plot,as.data.frame(pred_L$summary[,c(7,11,12)]))
             }
             names(predVals) <- c("age_plot","fit","lower","upper")
             plot(Lt ~ t, type = "n", ylab = "L(t)", xlab = "t(age)", xlim=c(min(age_plot),max(age_plot)),
                  main = "Non-linear least squares method")
             polygon(c(age_plot,rev(age_plot)),
                     c(predVals$lower,rev(predVals$upper)),
                     border = FALSE, col = adjustcolor("dodgerblue", alpha.f = 0.4))
             points(t, Lt)
             lines(age_plot, predict(nls_mod,newdata = data.frame(t=age_plot)))
           }
         },

         stop(paste("\n",method, "not recognised, possible options are \n",
                    "\"GullandHolt\", \"FordWalford\", \"Chapman\" \n",
                    "\"BertalanffyPlot\", and \"LSM\""))
         )

  res2 <- res
  if(exists("delta_t", where = tmp)) res2$delta_t <- delta_t
  if(exists("x", where = tmp)) res2$x <- x
  if(exists("y", where = tmp)) res2$y <- y
  if(exists("R2", where = tmp)) res2$r2 <- R2

  ret <- c(res2,list(mod = mod,
                    Linf = Linf,
                    K = K))

  if(exists("t0", where = tmp)) ret$t0 <- t0
  if(exists("nls_res", where = tmp)) ret$estimates <- nls_res

  return(ret)
}
