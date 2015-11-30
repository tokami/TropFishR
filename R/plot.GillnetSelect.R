#' @title Gillnet selectivity plot
#
#' @description  This function plots the selectivity estimates of the
#'   function \code{\link{GillnetSelect}}. Second sentence.
#'
#' @param x A list of the class \code{"GillnetSelect"} containing the results of the gillnet selectivity function.
#' @param ... normal parameters from plot function
#'
#' @examples
#' # load data
#' data(tilapia)
#'
#' # run model
#' output <- select(tilapia)
#'
#' # plot results
#' #plot(output)
#'
#' @details A function to plot the results of the gillnet selectivity estimation.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.GillnetSelect <- function(x,...){
  res <- x

  classes.num <- res$classes.num
  numNet1 <- res$numNet1
  numNet2 <- res$numNet2
  S1 <- res$S1
  S2 <- res$S2
  msNet1 <- res$msNet1
  msNet2 <- res$msNet2
  LmNet1 <- res$LmNet1
  LmNet2 <- res$LmNet2
  s2 <- res$stand.dev
  SNet1 <- res$SNet1
  SNet2 <- res$SNet2

  if(max(numNet2,na.rm=T) > max(numNet1,na.rm=T)){
    XX <- numNet2
  }else  XX <- numNet1

  classes.num.plot <- classes.num - 0.5
  xL <- seq(min(classes.num,na.rm=TRUE),max(classes.num,na.rm=TRUE),0.1)

  xlim <- c(min(xL,na.rm=TRUE),max(xL,na.rm=TRUE))
  ylim <- c(0, max(XX,na.rm=TRUE))
  plot.new()
  plot.window(xlim, ylim)

  #create plot
  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(classes.num.plot,XX,type='n', bty='n',xaxt ='n',
       ylab="numbers caught",xlab="fish length [cm]")
  lines(classes.num.plot,numNet1,type='s')
  axis(side=1,at=classes.num)
  lines(classes.num.plot,numNet2,type='s', lty=2)
  par(new=TRUE)
  plot(xL, exp(- ((xL - LmNet1)^2 / (2 * s2))), type = "l", col='darkgreen',lwd=2.5,
       axes=F,bty = "n", xlab = "", ylab = "")
  lines(xL, exp(- ((xL - LmNet2)^2 / (2 * s2))), type = "l", col='orange',lwd=2.5,
        bty = "n", xlab = "", ylab = "")
  axis(side=4, at = pretty(range(SNet1)))
  mtext("fractions retained", side=4, line=3)
  text(labels = paste("ms = ",msNet1,"cm"), x = LmNet1, y =
         ((exp(- ((LmNet1 - LmNet1)^2 / (2 * s2))))/1.1),
       col = 'darkgreen')
  text(labels = paste("ms = ",msNet2,"cm"), x = LmNet2, y =
         ((exp(- ((LmNet2 - LmNet2)^2 / (2 * s2))))/1.1),
       col = 'orange')
}
