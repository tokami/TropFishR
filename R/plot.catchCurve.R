#' @title Plotting catch curve
#
#' @description  This function plots the results from the \code{\link{catchCurve}} model.
#'
#' @param x A list of the class \code{"catchCurve"} containing the results of the
#'      catchCurve model.
#' @param plot.selec logical; if TRUE the regression line is plotted for not fully
#'      exploited length groups and the probability of capture is plotted
#' @param ... normal parameters from plot function
#'
#' @examples
#' \donttest{
#' #_______________________________________________
#' # Variable paramter system (with catch vector)
#' # based on length frequency data
#' # load data
#' data(goatfish)
#'
#' # run model
#' output <- catchCurve(goatfish, calc_ogive = TRUE)
#'
#' # plot results
#' plot(output)
#'
#' # based on age composition data
#' # load data
#' data(whiting)
#'
#' # run model
#' catchCurve(whiting, catch_column = 1)
#'
#' #_______________________________________________
#' # Catch Curve with estimation of selection ogive
#' # load data
#' data(synLFQ3)
#'
#' # run model
#' output <- catchCurve(synLFQ3, calc_ogive = TRUE)
#'
#' # plot results
#' plot(output,plot.selec = TRUE)
#' }
#'
#' @details A function to plot the results of the catchCurve model.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.catchCurve <- function(x, plot.selec = FALSE,...){
  pes <- x

  xlabel <- "Age [yrs]"
  if("t_midL" %in% names(pes)){
    xplot <- pes$t_midL
    xlabel <- "Relative age [yrs]"
    }else if("tplusdt_2" %in% names(pes)){
      xplot <- pes$tplusdt_2
    }else if("ln_Linf_L" %in% names(pes)){
        xplot <- pes$ln_Linf_L
        xlabel <- "ln(Linf - L)"
      }else if("classes.num" %in% names(pes)) xplot <- pes$classes.num

  if("lnC_dt" %in% names(pes)){
    yplot <- pes$lnC_dt
    ylabel <- "ln(C / dt)"
    }else if("lnC" %in% names(pes)){
      yplot <- pes$lnC
      ylabel <- "ln(C)"
    }else if("ln_C" %in% names(pes)){
        yplot <- pes$ln_C
        ylabel <- "ln C(L,Linf)"
    }
  if("ln_C" %in% names(pes) & "tplusdt_2" %in% names(pes)) ylabel <- "ln C(t, inf)"

  Z_lm1 <- pes$Z
  SE_Z_lm1 <- pes$se

  reg_int <- pes$reg_int


  #for plot
  minyplot <- ifelse(min(yplot,na.rm=TRUE) < 0, min(yplot,na.rm=TRUE),0)
  maxyplot <- max(yplot,na.rm=TRUE) + 1

  if(plot.selec){
    maxyplot <- ceiling(pes$intercept)

    op <- par(mfrow=c(2,1), xpd = FALSE,
              mar = c(1.2, 4, 1, 1) + 0.1,
              oma = c(6, 0.5, 1, 2) + 0.1)
    #final plot
    plot(x = xplot, y = yplot, ylim = c(minyplot,maxyplot),
         xlab = '', xaxt = 'n', ylab = ylabel,
         cex = 1.5)
    #par(new=T)
    points(x = xplot[reg_int[1]:reg_int[2]], y = yplot[reg_int[1]:reg_int[2]],
           pch = 19, col = 'blue', cex = 1.5)
    segments(x0=xplot[reg_int[1]],
             y0=yplot[reg_int[1]],
             x1=xplot[reg_int[2]],
             y1=yplot[reg_int[2]],
             col="blue",lwd = 1.7)
    segments(x0=0,
             y0=pes$intercept,
             x1=xplot[reg_int[1]],
             y1=yplot[reg_int[1]],
             col="blue",lwd = 1.7, lty = 2)
    mtext(side = 3, line = 0.3,  text = paste("Z =",round(Z_lm1,2),"+/-",
                                 round(SE_Z_lm1,2)), col = 'blue')

    plot(pes$Sest ~ xplot, type ='o', xlab = xlabel,
         ylab = "Probability of capture")
    points(y = 0.5, x=pes$t50,col='red',pch=16)
    segments(x0 = pes$t50, y0 = 0.5, x1 = pes$t50, y1 = 0, col='red',lty=2)
    op1 <- par(xpd=TRUE)
    text(y=-0.05, x=pes$t50, labels = "t50%", col = 'red')
    title(xlab = xlabel, outer = TRUE, line = 2)
    par(op)
    par(op1)


  }else {
    par(mfrow = c(1,1), mar = c(7, 5, 4, 5) + 0.3)
    #final plot
    plot(x = xplot, y = yplot, ylim = c(minyplot,maxyplot),
         xlab = xlabel, ylab = ylabel,
         cex = 1.5)
    par(new=T)
    points(x = xplot[reg_int[1]:reg_int[2]], y = yplot[reg_int[1]:reg_int[2]],
           pch = 19, col = 'blue', cex = 1.5)
    segments(x0=xplot[reg_int[1]],
             y0=yplot[reg_int[1]],
             x1=xplot[reg_int[2]],
             y1=yplot[reg_int[2]],
             col="blue",lwd = 1.7)
    mtext(side = 3, text = paste("Z =",round(Z_lm1,2),"+/-",
                                 round(SE_Z_lm1,2)), col = 'blue')
  }
}

