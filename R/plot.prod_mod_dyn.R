#' @title Plotting dynamic production models
#'
#' @description This function plots objects of the class "prod_mod_dyn".
#'
#' @param x a object of the class 'prod_mod_dyn',
#' @param ... optional parameters of plot function
#'
#' @examples
#' # load data
#' data(emperor)
#'
#' # run model
#' output <-  prod_mod_dyn(data = emperor)
#'
#' # plot output
#' plot(output)
#'
#'
#' @references
#' Dharmendra, D., Sólmundsson, J., 2005. Stock assessment of the offshore Mauritian banks
#' using dynamic biomass models and analysis of length frequency of the Sky Emperor
#' (\emph{Lethrinus mahsena}). Fisheries Training Program The United Nations University, 61
#'
#' @export

plot.prod_mod_dyn <- function(x,...){
  pes <- x
  Y <- pes$Y
  yrs <- pes$year
  Bvec <- pes$Bvec
  CPUE <- pes$CPUE
  Ihat <- pes$Ihat
  K <- pes$K
  r <- pes$r
  q <- pes$q
  MSY <- pes$MSY
  Bmsy <- pes$Bmsy
  Emsy <- pes$Emsy

  # plotting
  op <- par(mfrow=c(3,3))
  plot(yrs, Y, type="b", xlab="year", ylab="yield",
       main="Yield Trajectory NAZARETH",
       ylim=c(0, max(Y)*1.05))
  plot(yrs, Bvec, type="b", xlab="year", ylab="biomass",
       main="Biomass Trajectory", ylim=c(0, max(Bvec)*1.05))
  plot(Y, Bvec, xlab="yield", ylab="predicted Biomass",
       main="Corr. between yield and biomass")
  #cor(Bvec,Y)
  plot(Bvec, CPUE, xlab="Bvec-biomass", ylab="I-cpue",
       main="Corr. between biomass and CPUE")
  lines(lowess(Bvec, CPUE), col=2)
  plot(log(Bvec), log(CPUE), xlab="Bvec-biomass", ylab="I-cpue",
       main="Corr. between biomass and CPUE")
  plot(yrs,CPUE, type="b",xlab="Year", ylab="CPUE", ylim=c(0, max(CPUE)*1.05))
  plot(CPUE, Ihat, xlab="observed CPUE", ylab="Ihat",
       xlim=c(0, max(c(CPUE,Ihat))), ylim=c(0, max(c(CPUE,Ihat))))
  plot(yrs, CPUE, xlab="years", ylab="CPUE", ylim=c(0, max(CPUE)*1.05), type="b")
  lines(yrs, Ihat, col=2, type="b")
  #cor(Bvec, CPUE)

  ## the equilibrium yield
  Blevels <- seq(0,ceiling(K),10)
  EYlevels <- r*Blevels*(1-Blevels/K)
  EYlevels <- ifelse(EYlevels<0, 0, EYlevels)
  plot(Blevels, EYlevels, type="l" , xlab="Biomass",
       ylab="Yield", main="Yield Curve NAZARETH",
       xlim = c(0,max(c(Blevels,Bvec))),
       ylim = c(0,max(c(EYlevels,Y))))
  lines(Bvec,Y,type="b", col=2)
  abline(h=MSY,col=2)
  abline(v=Bmsy, col=2)


#   # this plot is not working properly
#   Elevels <- Blevels
#   #Elevels <- seq(0,10000,1000)
#   plot(Elevels, EYlevels, type="l")
#   abline(h=MSY,col=2)
#   abline(v=Emsy, col=2)

  par(op)

}




