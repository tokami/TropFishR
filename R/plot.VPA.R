#' @title VPA plot
#
#' @description  This function plots the survivors, catches, natural losses, and fishing
#'    mortality resulting from the \link{VPA} model.
#'
#' @param x list of the class \code{"VPA"} containing the results of the VPA model.
#' @param display_last_class logical; should last age/length class be displayed in graph?
#' @param ... standard parameters of \code{\link{barplot}}
#'
#' @examples
#' data(whiting)
#' output <- VPA(whiting, terminalF = 0.5, analysis_type = "VPA")
#' plot(output, display_last_class = FALSE)
#'
#' @details A function to plot the results of the VPA model.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.VPA <- function(x, display_last_class = TRUE, ...){
  pes <- x
  survivors <- pes$survivors
  FM_calc <- pes$FM_calc
  classes.num <- pes$classes.num
  df.VPAnew <- pes$plot_mat

  if("age" %in% names(pes)) xlabel <- "Age"
  if("midLengths" %in% names(pes)) xlabel <- "Midlength [cm]"

  #save x axis positions
  max_sur <- round(max(survivors,na.rm=TRUE),digits=0)
  dim_sur <- 10 ^ (nchar(max_sur)-1)
  max_FM <- ceiling(max(FM_calc,na.rm=TRUE))
  max_clas <- max(classes.num,na.rm=TRUE)
  par(new = FALSE)
  mids <- barplot(df.VPAnew, xlab="", ann=TRUE, plot = FALSE,
                  ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur))

  #create VPA plot
  dev.new()
  par(mar = c(7, 5, 4, 5) + 0.3)
  barplot(df.VPAnew,col=c('darkgreen','darkmagenta','gold2'),
          xlab = xlabel, ylab = "Population",xlim=c(0,ceiling(max(mids))),
          ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur),...)
  legend(x = mids[(which(classes.num == max_clas)-(length(classes.num)/4))],
         y = ceiling(max_sur/dim_sur)*dim_sur,
         legend = c(rownames(df.VPAnew),"fishing mortality"),
         col = c('darkgreen','darkmagenta','gold2','red'),xpd = TRUE,
         pch=c(rep(15,3),NA), lty = c(NA,NA,NA,1), lwd=2, seg.len = 0.9,
         pt.cex = 2, x.intersp = c(0.7,0.7,0.7,0.7), merge=TRUE,
         y.intersp = 1.2, box.lty=0, cex=0.8, xjust = -0.3, yjust = 0.7)
  par(new = TRUE, mar=c(7, 5, 4, 5) + 0.3)
  plot(mids, FM_calc, col='red',xlim=c(0,ceiling(max(mids))),
       type = "n",axes = FALSE, bty = "n", xlab = "", ylab = "", ann=TRUE)
  lines(x=mids,y=FM_calc,col='red',lwd=2)
  usr <- par("usr")
  par(usr=c(usr[1:2], 0, max_FM), mar=c(7, 5, 4, 5) + 0.3)
  axis(4, at = pretty(c(0,max_FM)))
  mtext("Fishing mortatlity", side=4, line=3)
}

