#' @title Plotting production models
#'
#' @description This function plots CPUE and yield values against fishing effort
#'      resulting from the production models (\link{prod_mod}).
#'
#' @param x a object of the class \code{"prod_mod"}
#' @param display_MSY logical; should MSY be displayed in the graph?
#' @param ... optional parameters of plot function
#'
#' @examples
#' data(trawl_fishery_Java)
#' output <-  prod_mod(data = trawl_fishery_Java)
#' plot(output, display_years = TRUE)
#'
#' @importFrom graphics abline layout legend lines par plot points segments text title
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
#' @method plot prod_mod
#' @export

plot.prod_mod <- function(x, display_MSY = TRUE, ...){
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
  CPUE.S <- pes$CPUE
  CPUE.F <- pes$ln_CPUE

  #Plot
  op <- par(mfrow=c(2,1), xpd = FALSE,
            mar = c(1.2, 4, 1.5, 1) + 0.1,
            oma = c(5, 0.5, 2, 2) + 0.1)
  layout(matrix(c(1,2), nrow = 2, byrow=TRUE), heights=c(4,4))

  # CPUE plot
  x <- seq(from = min(f),to = max(f),length.out = 500)
  y <- exp(a.F + b.F*x)
  plot(CPUE.S ~ f, xlab='', ylab='CPUE')
  lines(x = x, y = (a.S + x * b.S))
  #abline(a = a.S, b=b.S)
  lines(x = x, y = y, col = 'blue', lty = 5)

  # biomass plot
  x = seq(0,abs(a.S / b.S),length.out = 500)
  y.S = a.S*x + b.S*x^2
  y.F = x * exp(a.F + b.F*x)
  plot(x,rep(max(Y),length(x)),type='n',
       ylim=c(0,max(Y)*1.1), xlab = "", ylab = "Yield")
  points(Y ~ f)
  lines(x, y = y.S)
  lines(x, y = y.F,col='blue', lty = 5)
  if(display_MSY){
    #segments(x0=fMSY.S,x1=fMSY.S,y0=0,y1=MSY.S,lty=2)
    #segments(x0=-500,x1=fMSY.S,y0=MSY.S,y1=MSY.S,lty=2)
    points(fMSY.S,MSY.S,col='red', pch=18)
    text(fMSY.S, MSY.S,col='red',labels = expression(paste("MSY"[Schaefer])),
         adj = c(-0.05,1.1))
    #segments(x0=fMSY.F,x1=fMSY.F,y0=0,y1=MSY.F,lty=2,col='blue')
    #segments(x0=-500,x1=fMSY.F,y0=MSY.F,y1=MSY.F,lty=2,col='blue')
    points(fMSY.F,MSY.F, col='red', pch=18)
    text(fMSY.F, MSY.F, col='red',labels = expression(paste("MSY"[Fox])),
         adj = c(-0.05,1.1))
  }
  title(xlab = "Fishing effort", outer = TRUE, line = 2)
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(x="top", col=c('black','blue'), bty ='n', horiz = TRUE,
         xpd = TRUE,legend = c('Schaefer','Fox'), lty=c(1,5))
  par(op)
}
