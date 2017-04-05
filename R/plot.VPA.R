#' @title VPA plot
#'
#' @description  This function plots the survivors, catches, natural losses, and fishing
#'    mortality resulting from the \link{VPA} model.
#'
#' @param x list of the class \code{"VPA"} containing the results of the VPA model.
#' @param yaxis indicating which variable should be displayed on the y axis, either "numbers" or "biomass".
#' @param display_last_class logical; should last age/length class be displayed in graph?
#' @param xlabel Label of the x axis
#' @param ylabel1 Label of the first y axis
#' @param ylabel2 Label of the second y axis
#' @param ylim limits of y axis
#' @param ylim_FM limits of y axis of fishing mortality plot
#' @param plot.bars logical; should the barplot of survivors, nat.losses and catch
#'     be displayed? (Default: TRUE)
#' @param plot.FM logical; should the fishing mortality be displayed in the graph? (Default: TRUE)
#' @param plot.legend logical; should a legend be displayed in the graph? (Default: TRUE)
#' @param ... standard parameters of \code{\link{barplot}}
#'
#' @examples
#' data(whiting)
#' output <- VPA(whiting, terminalF = 0.5)
#' plot(output, display_last_class = FALSE)
#'
#' data(hake)
#' output <- VPA(hake, terminalE = 0.5)
#' plot_mat <- output$plot_mat[,-c(1,2)]  # remove first two length classes
#' class(plot_mat) <- "VPA"
#' plot(plot_mat, xlabel = "Midlengths [cm]")
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics axis barplot legend lines par plot mtext
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.VPA <- function(x,
                     yaxis = "numbers",
                     display_last_class = TRUE,
                     xlabel = NA,
                     ylabel1 = "Population",
                     ylabel2 = "Fishing mortality",
                     ylim = NA,
                     ylim_FM = NA,
                     plot.bars = TRUE,
                     plot.FM = TRUE,
                     plot.legend = TRUE,
                     ...){
  pes <- x
  if(!is.list(pes)){pes <- list(plot_mat = x)}
  survivors <- pes$plot_mat[1,]
  natLoss <- pes$plot_mat[2,]
  catch <- pes$plot_mat[3,]
  FM_calc <- pes$plot_mat[4,]
  if(dim(pes$plot_mat)[1] == 5){
    meanBodyWeight <- pes$plot_mat[5,]
  }
  if(dim(pes$plot_mat)[1] == 6){
    meanBodyWeight <- pes$plot_mat[5,]
    meanBiomassTon <- pes$plot_mat[6,]
  }
  classes.num <- as.numeric(colnames(pes$plot_mat))

  #put together in dataframe
  if(yaxis == "numbers"){
    df.VPAnew <- data.frame(survivors = survivors,
                            nat.losses = natLoss,
                            catch = catch)
  }
  if(yaxis == "biomass"){
    df.VPAnew <- data.frame(survivors = survivors * meanBodyWeight, #c(meanBiomassTon[-1],0) * 1000,
                            nat.losses = natLoss * meanBodyWeight,
                            catch = catch * meanBodyWeight)
  }

  #transpose matrix for barplot function
  df.VPAnew <- t(as.matrix(df.VPAnew))
  colnames(df.VPAnew) <- classes.num

  if(is.na(xlabel)){
    if("age" %in% names(pes)){
      xlabel <- "Age"
    }else if("midLengths" %in% names(pes)){
      xlabel <- "Midlength [cm]"
    }else{
      xlabel = "NA"
    }
  }else{
    xlabel = xlabel
  }


  #save x axis positions
  max_sur <- round(max(colSums(df.VPAnew),na.rm=TRUE),digits=0)
  if(is.na(ylim[1 ])){
    ylim = c(0, max_sur)
  }else{
    ylim  = ylim
  }
  dim_sur <- 10 ^ (nchar(max_sur)-1)
  max_FM <- ceiling(max(FM_calc,na.rm=TRUE))
  if(is.na(ylim_FM[1])){ylim_FM <- c(0,max_FM)}
  max_clas <- max(classes.num,na.rm=TRUE)
  par(new = FALSE)
  mids <- barplot(df.VPAnew, xlab="", ann=TRUE, plot = FALSE,
                  ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur))

  #create VPA plot
  #dev.new()
  if(par("mfrow")[1] == 1 & par("mfrow")[2] == 1){
    op <- par(mar = c(7, 5, 4, 5))
  }
  if(plot.bars){
    barplot(df.VPAnew,col=c('darkgreen','darkmagenta','gold2'),
            xlab="",
            ylab = ylabel1, xlim=c(0,ceiling(max(mids))),
            ylim = ylim, yaxs="i", ...)
    mtext(text = xlabel, side = 1, line = 2.5)
  }
  if(plot.legend){
    legend("topright",
           legend = c("Catch","Natural losses","Survivors","Fishing mortality"),
           col = c('gold2','darkmagenta','darkgreen','red'),xpd = TRUE,
           pch=c(rep(15,3),NA), lty = c(NA,NA,NA,1), lwd=2, seg.len = 0.9,
           pt.cex = 2, x.intersp = c(0.7,0.7,0.7,0.7), merge=TRUE,
           y.intersp = 1.2, box.lty=0, cex=0.8, xjust = -0.3, yjust = 0.7)
  }
  if(plot.FM){
    par(new = TRUE)
    plot(mids, FM_calc, col='red',xlim=c(0,ceiling(max(mids, na.rm = TRUE))),
         ylim=ylim_FM,
         type = "n", yaxs="i", axes = FALSE, bty = "n", xlab = "", ylab = "")
    lines(x=mids,y=FM_calc,col='red',lwd=2)
    axis(4, at = pretty(c(0,max_FM)),line = 1)
    mtext(ylabel2, side=4, line=3.5)
  }
  if(par("mfrow")[1] == 1 & par("mfrow")[2] == 1){
    par(op)
  }
}

