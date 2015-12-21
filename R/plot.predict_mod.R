#' @title Plotting prediction models yield per recruit and Thompson & Bell
#'
#' @description This function plots objects of the class "predict_mod", which are results
#'    of the function \code{\link{predict_mod}}.
#'
#' @param x a object of the class 'predict_mod',
#' @param ... optional parameters of plot function
#'
#' @examples
#'
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.predict_mod <- function(x,...){
  pes <- x

  #THOMPBELL
  if("totals" %in% names(pes)){


    #save x axis positions
    max_val <- round(max(pes$tot.V,na.rm=TRUE),digits=0)
    dim_val <- 10 ^ (nchar(max_val)-1)
    max_yiel <- round(max(pes$tot.Y,na.rm=TRUE),digits=0)
    dim_yiel <- 10 ^ (nchar(max_yiel)-1)
    max_bio <- round(max(pes$meanB,na.rm=TRUE),digits=0)
    dim_bio <- 10 ^ (nchar(max_bio)-1)

    #   max_val <- round(max(pes$tot.V,na.rm=TRUE),digits=0)
    #   dim_val <- 10 ^ (nchar(max_val)-1)

    op <- par(oma = c(1, 1, 1.5, 1),new=FALSE,mar = c(5, 4, 4, 6) + 0.3)
    plot(pes$Xfact,pes$tot.V, type ='o',ylab='Value',xlab='F-factor X',
         col ='darkorange', ylim = c(0,ceiling(max_val/dim_val)*dim_val),
         lwd=1.6)
    par(new=TRUE)
    plot(pes$Xfact,pes$tot.Y,type ='o',ylab='',xlab='',
         col='dodgerblue',lwd=1.6,axes=FALSE,
         ylim = c(0,ceiling(max_yiel/dim_yiel)*dim_yiel))
    axis(4,at=pretty(c(0,pes$tot.Y)))
    mtext("Yield", side=4, line=2)
    par(new=TRUE)
    plot(pes$Xfact,pes$meanB,type='o',axes=FALSE,ylab='',xlab='',
         col = 'darkgreen',lwd=1.6,
         ylim = c(0,ceiling(max_bio/dim_bio)*dim_bio))    # draw lines with small intervals: seq(0,max(),0.05) but y as to be dependent of x (formula of calculaiton of y)
    axis(4,at=pretty(c(0,pes$meanB)),line = 3)
    mtext("Biomass", side=4, line=5)

    par(oma = c(0, 0, 0, 0), new = TRUE)
    legend("top", c("value", "yield", "biomass"), xpd = TRUE,
           horiz = TRUE, inset = c(0, -0.1), bty = "n",lty = 1,seg.len = 0.7,
           col = c('darkorange','dodgerblue','darkgreen'), cex = 0.8,lwd=2,
           text.width=0.3,x.intersp=0.3)
    par(op)


  }

  # THOMP BELL WITH LC change
  if("mat_FM_Lc_com.C" %in% names(pes)){
    FM_change <- pes$FM_change
    Lc_change <- pes$Lc_tc_change
    mat_FM_Lc_com.Y <- pes$mat_FM_Lc_com.Y

    # colours for plot
    pal <- colorRampPalette(c(
      rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1)
    ))

    #plot
    image(x = FM_change,
          y = Lc_change,
          z = mat_FM_Lc_com.Y, col=pal(100),
          xlab = 'Fishing mortality', ylab = 'Lc')
    contour(x = FM_change,
            y = Lc_change,
            z = mat_FM_Lc_com.Y, add=TRUE)
    #mtext("Yield", line=0.5, side=3)

  }

  # YPR
  if("list_Lc_runs" %in% names(pes) |
     "list_tc_runs" %in% names(pes)){
    FM <- pes$FM
    if("tc" %in% names(pes)) tc_Lc <- pes$tc
    if("Lc" %in% names(pes)) tc_Lc <- pes$Lc
    if("list_tc_runs" %in% names(pes)) list_tc_Lc_runs <- pes$list_tc_runs
    if("list_Lc_runs" %in% names(pes)) list_tc_Lc_runs <- pes$list_Lc_runs

    #Plot
    op <- par(mfrow=c(1,1),new=F, mar = c(5, 4, 4, 4) + 0.3)

    #plot Y/R & B/R
    label <- ifelse("tc" %in% names(pes), "tc", "Lc")
    tc_Lc_start <- which(tc_Lc == max(tc_Lc,na.rm=T))
    offset_text <- list_tc_Lc_runs[[tc_Lc_start]]$Y_R[length(list_tc_Lc_runs[[tc_Lc_start]]$Y_R)] *
      0.02
    plot(list_tc_Lc_runs[[tc_Lc_start]]$FM, list_tc_Lc_runs[[tc_Lc_start]]$Y_R, type = 'l',
         ylab = "Y/R", xlab = "F",lty=tc_Lc_start,...)
    text(x = list_tc_Lc_runs[[tc_Lc_start]]$FM[length(list_tc_Lc_runs[[tc_Lc_start]]$FM)],
         y = (list_tc_Lc_runs[[tc_Lc_start]]$Y_R[length(list_tc_Lc_runs[[tc_Lc_start]]$Y_R)] +
                offset_text),
         labels = bquote(.(label)[.(names(list_tc_Lc_runs)[tc_Lc_start])]))
    seq_tc_Lc <- 1:length(list_tc_Lc_runs)
    seq_tc_Lc <- seq_tc_Lc[-tc_Lc_start]
    for(j in seq_tc_Lc){
      lines(list_tc_Lc_runs[[j]]$FM, list_tc_Lc_runs[[j]]$Y_R, type = 'l',
            ylab = "Y/R", xlab = "F",lty = j)
      text(x = list_tc_Lc_runs[[j]]$FM[length(list_tc_Lc_runs[[j]]$FM)],
           y = (list_tc_Lc_runs[[j]]$Y_R[length(list_tc_Lc_runs[[j]]$Y_R)] +
                  offset_text),
           labels = bquote(.(label)[.(names(list_tc_Lc_runs)[j])]))
    } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area
    par(new=T)
    plot(list_tc_Lc_runs[[1]]$FM, list_tc_Lc_runs[[1]]$B_R, type = 'l',
         axes = F, ylab = '', xlab ='', lty = tc_Lc_start,
         col = 'blue')
    axis(side = 4, at = pretty(range(list_tc_Lc_runs[[1]]$B_R)), col = 'blue',
         col.axis = 'blue')
    mtext(side = 4, text = "B/R", line = 3, col = 'blue')
    for(j in seq_tc_Lc){
      lines(list_tc_Lc_runs[[j]]$FM, list_tc_Lc_runs[[j]]$B_R, type = 'l',
            ylab = "Y/R", xlab = "F", col='blue', lty = j)
    } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area
    par(op)
  }
}


