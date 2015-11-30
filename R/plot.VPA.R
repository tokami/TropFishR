#' @title VPA plot
#
#' @description  This function plots the survivors, catches, natural losses, and fishing
#'    mortality resulting from the \link{VPA} model.
#'
#' @param x A list of the class \code{"VPA"} containing the results of the VPA model.
#' @param ... normal parameters from plot function
#'
#' @examples
#' data(whiting)
#' output <- VPA(whiting, terminalF = 0.5, analysis.type = "VPA")
#' plot(output)
#'
#' @details A function to plot the results of the VPA model.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.VPA <- function(x,...){
  res <- x
  survivors <- res$survivors
  FM <- res$FM
  classes.num <- res$classes.num
  df.VPAnew <- res$plot_mat

  #save x axis positions
  max_sur <- round(max(survivors,na.rm=TRUE),digits=0)
  dim_sur <- 10 ^ (nchar(max_sur)-1)
  max_FM <- ceiling(max(FM,na.rm=TRUE))
  max_clas <- max(classes.num,na.rm=TRUE)
  par(new = FALSE)
  mids <- barplot(df.VPAnew, xlab="", ann=TRUE,
                  ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur))

  #create VPA plot
  par(mar = c(5, 4, 4, 4) + 0.3)
  barplot(df.VPAnew,col=c('darkgreen','purple','yellow'),
          xlab = "Age", ylab = "Population",xlim=c(0,ceiling(max(mids))),
          ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur), ...)
  legend(x=mids[(which(classes.num == max_clas)-2)],
         y = ceiling(max_sur/dim_sur)*dim_sur,
         legend = c(rownames(df.VPAnew),"Fishing mortality"),
         col = c('darkgreen','purple','yellow','red'),xpd = TRUE,
         pch=c(rep(15,3),NA), lty = c(NA,NA,NA,1), lwd=2,seg.len = 0.3,
         pt.cex = 2, x.intersp = c(0.3,0.3,0.3,0.3),merge=TRUE,
         y.intersp = 0.6, box.lty=0,cex=0.8,xjust = 0,yjust = 0.8)
  par(new = TRUE,mar = c(5, 4, 4, 4) + 0.3)
  plot(mids, FM, col='red',xlim=c(0,ceiling(max(mids))),
       type = "n",axes = FALSE, bty = "n", xlab = "", ylab = "",ann=TRUE)
  lines(x=mids,y=FM,col='red',lwd=2)
  usr <- par("usr")
  par(usr=c(usr[1:2], 0, max_FM))
  axis(4,at=pretty(c(0,max_FM)))
  mtext("fishing mortatlity", side=4, line=3)
}
