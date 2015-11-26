#' @title Trawl net selectivity plot
#
#' @description  This function estimates the selecitvity of trawl nets.
#'    First sentence. second sentence.
#'
#' @param x results of trawl select function
#' @param ... parameter of plot function
#'
#' @examples
#' data(data_TrawlSelect)
#' output <- TrawlSelect(data_TrawlSelect)
#' plot(output)
#'
#' @details To calculate selection factor (SF), L25, L50 and L75 for trawl nets /fisheries.
#'
#' @references xx
#'
#'
#' @export


plot.TrawlSelect <- function(x,...){
  res <- x

  classes.num <- res$classes.num
  SLobs <- res$SLobs
  S1 <- res$S1
  S2 <- res$S2
  L25 <- res$L25
  L50 <- res$L50
  L75 <- res$L75

  xL <- seq(min(classes.num,na.rm=TRUE),max(classes.num,na.rm=TRUE),0.1)

  xlim <- c(min(xL,na.rm=TRUE),max(xL,na.rm=TRUE))
  ylim <- c(0, max(SLobs,na.rm=TRUE))
  plot.new()
  plot.window(xlim, ylim)

  plot(SLobs ~ classes.num,
       xlab = "L [cm]", ylab = "fraction retained", pch = 16 , cex =1.4)
  lines(x = xL, y = 1/(1 + exp(S1 - S2 * xL)), col='blue', lwd = 1.5)
  segments(x0 = L25, x1 = L25, y0 = 0,
           y1 = 1/(1 + exp(S1 - S2 * L25)),col= 'gray40',lty = 2, lwd=1.5)
  segments(x0 = 0, x1 = L25, y0 = 1/(1 + exp(S1 - S2 * L25)),
           y1 = 1/(1 + exp(S1 - S2 * L25)),col= 'gray40',lty = 2, lwd=1.5)
  text(labels = "L25%", x = L25, y = (1/(1 + exp(S1 - S2 * L25)))/3,
       col = 'gray40')
  segments(x0 = L50, x1 = L50, y0 = 0,
           y1 = 1/(1 + exp(S1 - S2 * L50)),col= 'gray40',lty = 2, lwd=1.5)
  segments(x0 = 0, x1 = L50, y0 = 1/(1 + exp(S1 - S2 * L50)),
           y1 = 1/(1 + exp(S1 - S2 * L50)),col= 'gray40',lty = 2, lwd=1.5)
  text(labels = "L50%", x = L50, y = (1/(1 + exp(S1 - S2 * L50)))/3,
       col = 'gray40')
  segments(x0 = L75, x1 = L75, y0 = 0,
           y1 = 1/(1 + exp(S1 - S2 * L75)),col= 'gray40',lty = 2, lwd=1.5)
  segments(x0 = 0, x1 = L75, y0 = 1/(1 + exp(S1 - S2 * L75)),
           y1 = 1/(1 + exp(S1 - S2 * L75)),col= 'gray40',lty = 2, lwd=1.5)
  text(labels = "L75%", x = L75, y = (1/(1 + exp(S1 - S2 * L75)))/3,
       col = 'gray40')

}
