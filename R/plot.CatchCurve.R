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

plot.CatchCurve <- function(x,...){
  res <- x

  if("t_midL" %in% names(res)){
    xplot <- res$t_midL
    }else if("tplusdt_2" %in% names(res)){
      xplot <- res$tplusdt_2
      }else if("classes.num" %in% names(res)) xplot <- res$classes.num

  if("lnC_dt" %in% names(res)){
    yplot <- res$lnC_dt
    ylabel <- "ln(C / dt)"
    }else if("lnC" %in% names(res)){
      yplot <- res$lnC
      ylabel <- "ln(C)"
      }

  Z_lm1 <- res$Z
  SE_Z_lm1 <- res$se

  reg_int <- res$reg_int


  #for plot
  minyplot <- ifelse(min(yplot,na.rm=TRUE) < 0, min(yplot,na.rm=TRUE),0)
  maxyplot <- max(yplot,na.rm=TRUE) + 1

  #final plot
  plot(x = xplot, y = yplot, ylim = c(minyplot,maxyplot),
       xlab = "Relative age [yrs]", ylab = ylabel,
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
