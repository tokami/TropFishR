#' @title Plotting production models
#'
#' @description This function plots objects of the class "prod_mod".
#'
#' @param x a object of the class 'prod_mod',
#' @param ... optional parameters of plot function
#'
#' @examples
#' # load data
#' data(trawl_fishery_Java)
#'
#' # run model
#' output <-  prod_mod(data = trawl_fishery_Java)
#'
#' # plot output
#' plot(output)
#'
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
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.prod_mod <- function(x,...){
  pes <- x
  f <- pes$f
  Y <- pes$Y
  fMSY.S <- pes$Schaefer_fMSY
  MSY.S <- pes$Schaefer_MSY
  fMSY.F <- pes$Fox_fMSY
  MSY.F <- pes$Fox_MSY
  a.F <- pes$Fox_lm[1]
  b.F <- pes$Fox_lm[2]
  a.S <- pes$Schaefer_lm[1]
  b.S <- pes$Schaefer_lm[2]
  CPUE.S <- pes$Schaefer_CPUE

  #Plot
  op <- par(mfrow=c(2,1),new=F, mar=c(5, 4, 4, 2) + 0.1)
  x = seq(min(f),max(f),1)
  y = exp(a.F + b.F*x)
  plot(CPUE.S ~ f, xlab='Fishing effort', ylab='CPUE')
  abline(a = a.S, b=b.S)
  lines(x, y = y, col = 'blue', lty = 5)
  legend("topright", col=c('black','blue'), bty ='n',
         legend = c('Schaefer','Fox'), lty=c(1,5))


  x = seq(0,abs(a.S / b.S),1)
  y.S = a.S*x + b.S*x^2
  y.F = x * exp(a.F + b.F*x)
  par(new=F)
  plot(x,rep(max(Y),length(x)),type='n',
       ylim=c(0,max(Y)), xlab = "Fishing effort", ylab = "Yield")
  points(Y ~ f)
  lines(x, y = y.S)
  lines(x, y = y.F,col='blue', lty = 5)
  segments(x0=fMSY.S,x1=fMSY.S,y0=0,y1=MSY.S,lty=2)
  segments(x0=0,x1=fMSY.S,y0=MSY.S,y1=MSY.S,lty=2)
  points(fMSY.S,MSY.S,col='red',pch=16)
  segments(x0=fMSY.F,x1=fMSY.F,y0=0,y1=MSY.F,lty=2,col='blue')
  segments(x0=0,x1=fMSY.F,y0=MSY.F,y1=MSY.F,lty=2,col='blue')
  points(fMSY.F,MSY.F,col='red',pch=16)
  legend("bottomleft", col=c('black','blue'), bty ='n',
         legend = c('Schaefer','Fox'), lty=c(1,5))

  par(op)
}

