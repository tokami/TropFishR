#' @title Selectivity plot
#
#' @description  This function plots the selectivity estimates of the
#'   function \code{\link{select}}.
#'
#' @param x a list of the class \code{"select"} containing the results of the
#'   gillnet selectivity function.
#' @param regression_fit logical; indicating if a plot with the fit of the regression
#'    line should be displayed
#' @param ... additional parameters of the \link{plot} function
#'
#' @examples
#' data(tilapia)
#' output <- select(tilapia, plot = FALSE)
#' plot(output)
#'
#' data(bream)
#' output <- select(bream, plot = FALSE)
#' plot(output)
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.select <- function(x, regression_fit = FALSE, ...){
  res <- x
  classes.num <- res$classes.num

  if(res$type == "gillnet"){
    numNet1 <- res$CatchPerNet_mat[,1]
    numNet2 <- res$CatchPerNet_mat[,2]
    msNet1 <- res$meshSizes[1]
    msNet2 <- res$meshSizes[2]
    LmNet1 <- res$LmNet1
    LmNet2 <- res$LmNet2
    s2 <- res$stand.dev
    SNet1 <- res$SNet1
    SNet2 <- res$SNet2
    reg.coeffs <- res$reg.coeffs
    lnNet2_Net1 <- res$lnNet2_Net1

    if(max(numNet2,na.rm=T) > max(numNet1,na.rm=T)){
      XX <- numNet2
    }else  XX <- numNet1

    classes.num.plot <- classes.num - 0.5
    xL <- seq(min(classes.num,na.rm=TRUE),max(classes.num,na.rm=TRUE),0.1)

    dev.new()
    #create plot
    if(regression_fit){
      op <- par(mfrow = c(2,1), xpd = FALSE,
      mar = c(4, 4, 2, 3) + 0.1,
      oma = c(3, 0.5, 1, 2) + 0.1)

      plot(classes.num, lnNet2_Net1, xlab = "Length groups", ylab = "ln(Nnet2/Nnet1)",
           ylim = c(min(lnNet2_Net1,na.rm=TRUE)*1.05,
                    max(lnNet2_Net1,na.rm=TRUE)*1.05))
      abline(a = reg.coeffs[1], b = reg.coeffs[2])

    }else op <- par(mfrow = c(1,1), mar = c(5, 5, 3, 3),
                    oma = c(3, 0.5, 1, 2))

    plot(classes.num.plot,XX,type='n', bty='n',xaxt ='n',
         ylab="Numbers caught",xlab="Fish length")
    lines(classes.num.plot,numNet1,type='s')
    axis(side=1,at=classes.num)
    lines(classes.num.plot,numNet2,type='s', lty=2)
    par(new=TRUE)
    plot(xL, exp(- ((xL - LmNet1)^2 / (2 * s2))), type = "l", col='darkgreen',lwd=2.5,
         axes=F,bty = "n", xlab = "", ylab = "")
    lines(xL, exp(- ((xL - LmNet2)^2 / (2 * s2))), type = "l", col='orange',lwd=2.5,
          bty = "n", xlab = "", ylab = "")
    axis(side=4, at = pretty(range(SNet1)))
    mtext("Fractions retained", side=4, line=3)
    text(labels = paste("ms = ",msNet1,"cm"), x = LmNet1, y =
           ((exp(- ((LmNet1 - LmNet1)^2 / (2 * s2))))/0.98),
         col = 'darkgreen',xpd = TRUE)
    text(labels = paste("ms = ",msNet2,"cm"), x = LmNet2, y =
           ((exp(- ((LmNet2 - LmNet2)^2 / (2 * s2))))/0.98),
         col = 'orange',xpd = TRUE)
    par(op)
  }

  if(res$type == "trawl_net"){
    reg.coeffs <- res$reg.coeffs
    S1 <- res$S1
    S2 <- res$S2
    SLobs <- res$SLobs
    L25 <- res$L25
    L50 <- res$L50
    L75 <- res$L75
    lnSL <- res$lnSL

    xL <- seq(min(classes.num,na.rm=TRUE),max(classes.num,na.rm=TRUE),0.1)

    dev.new()
    if(regression_fit){
      op <- par(mfrow = c(2,1), xpd = FALSE,
                mar = c(4, 4, 2, 2) + 0.1,
                oma = c(3, 0.5, 1, 1) + 0.1)
      plot(classes.num, lnSL, xlab = "Length groups", ylab = "ln(Sl)",
           ylim = c(min(lnNet2_Net1,na.rm=TRUE)*1.05,
                    max(lnNet2_Net1,na.rm=TRUE)*1.05))
      abline(a=reg.coeffs[1], b=reg.coeffs[2])

    }else op <- par(mfrow = c(1,1), mar = c(5, 5, 3, 3),
                    oma = c(3, 0.5, 1, 2))

    plot(SLobs ~ classes.num,
         xlab = "Fish length", ylab = "Fraction retained", pch = 16 , cex =1.4)
    lines(x = xL, y = 1/(1 + exp(S1 - S2 * xL)), col='blue', lwd = 1.5)
    segments(x0 = L25, x1 = L25, y0 = -1,
             y1 = 1/(1 + exp(S1 - S2 * L25)),col= 'gray40',lty = 2, lwd=1.5)
    segments(x0 = 0, x1 = L25, y0 = 1/(1 + exp(S1 - S2 * L25)),
             y1 = 1/(1 + exp(S1 - S2 * L25)),col= 'gray40',lty = 2, lwd=1.5)
    text(labels = "L25%", x = L25, y = (1/(1 + exp(S1 - S2 * L25)))/3,
         col = 'gray40')
    segments(x0 = L50, x1 = L50, y0 = -1,
             y1 = 1/(1 + exp(S1 - S2 * L50)),col= 'gray40',lty = 2, lwd=1.5)
    segments(x0 = 0, x1 = L50, y0 = 1/(1 + exp(S1 - S2 * L50)),
             y1 = 1/(1 + exp(S1 - S2 * L50)),col= 'gray40',lty = 2, lwd=1.5)
    text(labels = "L50%", x = L50, y = (1/(1 + exp(S1 - S2 * L50)))/3,
         col = 'gray40')
    segments(x0 = L75, x1 = L75, y0 = -1,
             y1 = 1/(1 + exp(S1 - S2 * L75)),col= 'gray40',lty = 2, lwd=1.5)
    segments(x0 = 0, x1 = L75, y0 = 1/(1 + exp(S1 - S2 * L75)),
             y1 = 1/(1 + exp(S1 - S2 * L75)),col= 'gray40',lty = 2, lwd=1.5)
    text(labels = "L75%", x = L75, y = (1/(1 + exp(S1 - S2 * L75)))/3,
         col = 'gray40')
    par(op)
  }
}

