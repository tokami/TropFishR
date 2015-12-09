#' @title ELEFAN plot
#
#' @description  This function plots the growth curves calculated by means of ELEFAN
#'    (\code{\link{ELEFAN}}).
#'
#' @param x A list of the class \code{"ELEFAN"} containing the results of ELEFAN
#' @param ... Default parameter options from plot function
#'
#' @examples
#' \donttest{
#' #load data
#'   Linf.p = 45.2
#'   K.p = 0.6
#'   tshift.p = 360/365
#'
#' plot(output, Linf= , K=, )  # with entering tshift? or it looks up tshift from time matrix but what happens if user choose L and K values outside of time mat???
#'
#' }

#'
#' @details This function draws a selectivity plot for the object class
#'    \code{"ELEFAN"}, which is created by applying ELEFAN
#'    \code{\link{ELEFAN}}.
#'
#' @references
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Pauly, D., 1987. A review of the ELEFAN system for analysis of length-frequency data in
#' fish and aquatic invertebrates. \emph{ICLARM Conf. Proc.}, (13):7-34
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.ELEFAN <- function(x,...
                        ){
  res <- x

  # Plot with rearranged histogramms
  # axes
  xplot <- c(min(t.years.fr),(max(t.years.fr)))
  yplot <- c(classes[1]-interval,classes[length(classes)]+interval)

  par(mar = c(5, 5, 1, 1) + .1)
  plot(xplot, yplot, type = "n",axes = FALSE,
       ann = FALSE,  xaxs = "i", yaxs = "i")
  text(days.years, par("usr")[3] - 0.7,
       labels = dates.all, srt = 45, pos = 1, xpd = TRUE)
  axis(2, cex.axis = 1.2)
  mtext(side = 2, outer = F, line = 3.5, "Length [cm]", cex = 1.2)
  box( col = "gray40") #bty = "L"

  #Histograms
  # positive histograms or peaks
  y.b <- classes - interval / 2
  y.t <- classes + interval / 2
  for(histi in 1:dim(catch)[2]){
    x.l <- rep(days.years[histi],length(classes))
    x.x <- catch.aAF[,histi]
    x.x[which(x.x < 0)] <- 0
    x.r <- x.l - x.x * 0.05   ### make this number dependent on data!?!?
    rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
         col = "gray40", border = "black")
  }

  # negative histograms or troughs
  for(histi in 1:dim(catch)[2]){
    x.l <- rep(days.years[histi],length(classes))
    x.x <- catch.aAF[,histi]
    x.x[which(x.x > 0)] <- 0
    x.r <- x.l - x.x * 0.05   ### make this number dependent on data!?!?
    rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
         col = "gray80", border = "black")
  }



  # growth curves
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p)))    ### choose correct t -> tshift
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')
  #t.years.frA <- seq(xplot[1],xplot[2],1/365)
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p - 1)))
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p + 1)))
  lines(y = Lt, x = t.years.fr , lty=2, col='blue')
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p + 2)))
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')
  Lt <- Linf.p * (1 - exp(-K.p * (t.years.fr - tshift.p + 3)))
  lines(y = Lt, x = t.years.fr, lty=2, col='blue')



  #for plotting
  xplot <- c(min(t.years.fr),(max(t.years.fr)))     # display months on x axis
  yplot <- c(classes[1]-interval,classes[length(classes)]+interval)

  # Plot with rearranged histogramms
  par(mar = c(5, 5, 1, 1) + .1)
  plot(xplot, yplot, type = "n",axes = FALSE,
       ann = FALSE,  xaxs = "i", yaxs = "i")
  text(days.years, par("usr")[3] - 0.7,
       labels = dates.all, srt = 45, pos = 1, xpd = TRUE)
  axis(2, cex.axis = 1.2)
  mtext(side = 2, outer = F, line = 3.5, "Length [cm]", cex = 1.2)
  box( col = "gray40") #bty = "L"

  #Histograms
  # true distributions
  for(histi in 1:dim(catch)[2]){
    x.r <- rep(days.years[histi],length(classes))
    x.x <- catch[,histi]
    x.l <- x.r - x.x * 0.005   ### make this number dependent on data!?!?
    rect(xleft = x.l, ybottom = y.b, xright = x.r, ytop = y.t,
         col = "gray80", border = "black")
  }




  #   for(xyx in 1:length(Lt.list)){
  #     lines(y = Lt.list[[xyx]], x = t, lty=2, col='darkgreen')
  #   }

}




