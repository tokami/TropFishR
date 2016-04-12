#' @title Plotting time series production models
#'
#' @description This function plots objects of the class "prod_mod_ts".
#'
#' @param x a object of the class "prod_mod_ts",
#' @param correlation_plots logical; indicating if correlation plots should be displayed
#' @param ... additional parameters of the \link{plot} function
#'
#' @examples
#' data(emperor)
#' output <-  prod_mod_ts(emperor, method  = "Fox")
#' plot(output, correlation_plots = TRUE)
#'
#'
#' @references
#' Dharmendra, D., Solmundsson, J., 2005. Stock assessment of the offshore Mauritian banks
#' using dynamic biomass models and analysis of length frequency of the Sky Emperor
#' (\emph{Lethrinus mahsena}). Fisheries Training Program The United Nations University, 61
#'
#' @export

plot.prod_mod_ts <- function(x, correlation_plots = FALSE, ...){
  pes <- x
  Y <- pes$Y
  yrs <- pes$year
  Bvec <- pes$Bvec
  CPUE <- pes$CPUE
  CPUE_hat <- pes$CPUE_hat
  K <- pes$K
  r <- pes$r
  q <- pes$q
  MSY <- pes$MSY
  Bmsy <- pes$Bmsy
  Emsy <- pes$Emsy
  method <- pes$method

  # plotting
  dev.new()
  op <- par(mfrow=c(2,2), xpd = FALSE,
            mar = c(4, 4, 3, 1) + 0.1,
            oma = c(4, 0.5, 1, 2) + 0.1)
  layout(matrix(c(1,2,3,4), nrow = 2, byrow=TRUE))

  # Yiel trajectory
  plot(yrs, Y, type="b", xlab="", ylab="Yield",
       main="Yield Trajectory",
       ylim=c(0, max(Y, na.rm = TRUE)*1.05))

  # Biomass trajectory
  plot(yrs, Bvec, type="b", xlab="", ylab="Biomass",
       main="Biomass Trajectory", ylim=c(0, max(Bvec, na.rm = TRUE)*1.05))

  # CPUE trajectory
  plot(yrs, CPUE, xlab="", ylab="CPUE",
       ylim=c(0, max(CPUE, na.rm = TRUE)*1.05), type="b",
       main = "CPUE Trajectory")
  lines(yrs, CPUE_hat, col=2, type="b")
  legend(x="bottomright",legend = c("observed", "predicted"),lty = 1, pch = 1,xpd = TRUE,
         col = c(1,2), cex = 0.8, bty = 'n', y.intersp = 0.8, x.intersp = 0.5)


  ## the equilibrium yield
  Blevels <- seq(0,ceiling(K),10)
  if(method == "Schaefer") EYlevels <- r*Blevels*(1-Blevels/K)
  if(method == "Fox") EYlevels <- r*Blevels*log(K/Blevels)
  EYlevels <- ifelse(EYlevels<0, 0, EYlevels)
  plot(Blevels, EYlevels, type="l" , xlab="Biomass",
       ylab="Yield", main="Yield Curve",
       xlim = c(0,max(c(Blevels,Bvec), na.rm = TRUE)),
       ylim = c(0,max(c(EYlevels,Y), na.rm = TRUE)))
  lines(Bvec,Y,type="b", col=2)
  abline(h=MSY,col=2)
  abline(v=Bmsy, col=2)

  par(op)

  # Validating Fit
  if(correlation_plots){
    # plotting
    dev.new()
    op <- par(mfrow=c(2,2), xpd = FALSE,
              mar = c(4, 4, 4, 1) + 0.1,
              oma = c(4, 0.5, 1, 2) + 0.1)
    layout(matrix(c(1,2,3,4), nrow = 2, byrow=TRUE))

    # Correlation between yield and biomass
    plot(Y, Bvec, xlab="Yield", ylab="Biomass",
         main="Corr. between yield and biomass")
    lines(lowess(Y, Bvec), col=2)

    # Correlation between biomass and CPUE
    plot(Bvec, CPUE, xlab="Biomass", ylab="CPUE",
         main="Corr. between biomass and CPUE")
    lines(lowess(Bvec, CPUE), col=2)

    # Correlation between biomass and CPUE
    plot(log(Bvec), log(CPUE), xlab="Biomass", ylab="CPUE",
         main="Corr. between biomass and CPUE")
    lines(lowess(log(Bvec), log(CPUE)), col=2)

    # observed vs predicted CPUE
    plot(CPUE, CPUE_hat, xlab="observed CPUE", ylab="predicted CPUE",
         main="Corr. observed and predicted CPUE")
    lines(lowess(CPUE, CPUE_hat), col=2)

    par(op)
  }
}
