#' @title Catchcurve plot
#
#' @description  This function plots the results from the \link{CatchCurve} model.
#'
#' @param x A list of the class \code{"CatchCurve"} containing the results of the CatchCurve model.
#' @param ... normal parameters from plot function
#'
#' @examples
#' \donttest{
#' data(whiting)
#' output <- CatchCurve(whiting, catch_column = 1)
#' plot(output)
#' }
#'
#' @details A function to plot the results of the CatchCurve model.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export
x <- output
plot.CatchCurve <- function(x,...){
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
