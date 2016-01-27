#' @title Plotting prediction models yield per recruit and Thompson & Bell
#'
#' @description This function plots objects of the class "predict_mod", which are results
#'    of the function \code{\link{predict_mod}}.
#'
#' @param x a object of the class 'predict_mod'
#' @param type default = "ypr" or "Isopleth" for isopleth plot
#' @param xaxis1 default = "FM" which x axis should be plotted? FM or E?
#' @param yaxis1  default = "Y_R" for ypr plot, which first y axis should be plotted? Y_R or Y_R.rel
#' @param yaxis2 default = "B_R" for ypr plot, which second y axis should be plotted? B_R or B_R.percent
#' @param ... optional parameters of plot function
#'
#' @examples
#' # YPR
#' # age structured data
#' # Nemipterus marginatus
#' threadfin <- list(Winf = 286,K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4)
#'
#' # run model
#' predict_mod(threadfin, FM_change = seq(0,6,0.1),
#'    Lc_tc_change = seq(0.2,1,0.2), type = 'ypr')  #where it is maximal  = MSY
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.predict_mod <- function(x, type = 'ypr', xaxis1 = "FM",
                             yaxis1 = "Y_R", yaxis2 = "B_R",...){
  pes <- x <- output

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
         ylim = c(0,ceiling(max_yiel/dim_yiel) * dim_yiel))
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

  # THOMP BELL WITH LC change - ISOPLETHS
  if("mat_FM_Lc_com.C" %in% names(pes)){
    FM_change <- pes$FM_change
    Lc_change <- pes$Lc_tc_change
    mat_FM_Lc_com.Y <- pes$mat_FM_Lc_com.Y

    # colours for plot
    pal <- colorRampPalette(rev(c(
      rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1))))

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

    # necessary calculations
    FM <- pes$FM
    if("tc" %in% names(pes)) tc_Lc <- pes$tc
    if("Lc" %in% names(pes)) tc_Lc <- pes$Lc
    if("list_tc_runs" %in% names(pes)) list_tc_Lc_runs <- pes$list_tc_runs
    if("list_Lc_runs" %in% names(pes)) list_tc_Lc_runs <- pes$list_Lc_runs

    p.yield <- yaxis1
    p.FE <- xaxis1
    p.B <- yaxis2
    xlabel1 <- ifelse(xaxis1 == "FM", "Fishing mortality", "Exploitation rate")
    ylabel1 <- ifelse(yaxis1 == "Y_R", "Y/R", "rel. Y/R")
    ylabel2 <- ifelse(yaxis2 == "B_R", "B/R", "B/R [%]")

    df_Es <- pes$df_Es
    if(xaxis1 == "FM"){
      N01 <- df_Es$F01
      N05 <- df_Es$F05
      Nmax <- df_Es$Fmax
      legend.lab <- c("F0.1","F0.5","Fmax")
    }else{
      N01 <- df_Es$E01
      N05 <- df_Es$E05
      Nmax <- df_Es$Emax
      legend.lab <- c("E0.1","E0.5","Emax")
    }

    #standard plot (ypr vs. E or FM)
    if(type == 'ypr'){
      #Plot
      op <- par(mfrow=c(1,1),new=F, mar = c(5, 4, 4, 4) + 0.3)

      #plot Y/R & B/R
      label <- ifelse("tc" %in% names(pes), "tc", "Lc")
      tc_Lc_start <- which(tc_Lc == max(tc_Lc,na.rm=T))
      p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
      py <- p.dat[,which(names(p.dat) == p.yield)]
      px <- p.dat[,which(names(p.dat) == p.FE)]
      offset_text <- py[length(py)] * 0.02

      plot(px, py, type = 'l',
           ylab = ylabel1, xlab = xlabel1, lty=tc_Lc_start)  #,...
      text(x = px[length(px)],
           y = (py[length(py)] +
                  offset_text),
           labels = bquote(.(label)[.(names(list_tc_Lc_runs)[tc_Lc_start])]))
      seq_tc_Lc <- 1:length(list_tc_Lc_runs)
      seq_tc_Lc <- seq_tc_Lc[-tc_Lc_start]
      for(j in seq_tc_Lc){
        p.dat <- list_tc_Lc_runs[[j]]
        py <- p.dat[,which(names(p.dat) == p.yield)]
        px <- p.dat[,which(names(p.dat) == p.FE)]

        lines(px, py, type = 'l',
              ylab = ylabel1, xlab = xlabel1,lty = j)
        text(x = px[length(px)],
             y = (py[length(py)] +
                    offset_text),
             labels = bquote(.(label)[.(names(list_tc_Lc_runs)[j])]))
      } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area

      # reference points
      if(length(N01) == 1){
        p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
        py <- p.dat[,which(names(p.dat) == p.yield)]
        px <- p.dat[,which(names(p.dat) == p.FE)]
        # F or  E 0.1
        segments(x0 = -1, x1 = N01, y0 = py[which(px == N01)], y1 = py[which(px == N01)],
                 col= 'darkgreen',lty = 2, lwd=1.5)
        segments(x0 = N01, x1 = N01, y0 = -1, y1 = py[which(px == N01)],
                 col= 'darkgreen',lty = 2, lwd=1.5)
        # F or E max
        segments(x0 = -1, x1 = Nmax, y0 = py[which(px == Nmax)], y1 = py[which(px == Nmax)],
                 col= 'goldenrod1',lty = 2, lwd=1.5)
        segments(x0 = Nmax, x1 = Nmax, y0 = -1, y1 = py[which(px == Nmax)],
                 col= 'goldenrod1',lty = 2, lwd=1.5)
        # Legend
        legend("top", legend = legend.lab, xpd = TRUE, horiz = TRUE,
               inset = c(0,0), bty = "n", lty = 2, col = c("darkgreen","red","goldenrod1"),
               seg.len = 1,pt.cex = 2, x.intersp = c(0.7,0.7,0.7),merge=TRUE,
               y.intersp = -2, box.lty=0,cex=0.8, lwd =2)

        # current Exploitation rate or fishing mortality
        currents <- pes$currents
        if(!is.na(currents$curr.E)){
          px <- ifelse(p.FE == "FM",currents$curr.F, currents$curr.E)
          py <- ifelse(p.yield == "Y_R",currents$curr.YPR,currents$curr.YPR.rel)
          points(px,py, pch = 16)
        }
      }

      # Biomass per recruit
      par(new=T)
      px <- list_tc_Lc_runs[[1]][,which(names(list_tc_Lc_runs[[1]]) == p.FE)]
      py <- list_tc_Lc_runs[[1]][,which(names(list_tc_Lc_runs[[1]]) == p.B)]
      plot(px, py, type = 'l',
           axes = F, ylab = '', xlab ='', lty = tc_Lc_start,
           col = 'blue')
      axis(side = 4, at = pretty(range(py, na.rm= TRUE)), col = 'blue',
           col.axis = 'blue')
      mtext(side = 4, text = ylabel2, line = 3, col = 'blue')
      for(j in seq_tc_Lc){
        p.dat <- list_tc_Lc_runs[[j]]
        py <- p.dat[,which(names(p.dat) == p.B)]
        px <- p.dat[,which(names(p.dat) == p.FE)]

        lines(px, py, type = 'l',
              ylab = ylabel1, xlab = xlabel1, col='blue', lty = j)
      }
      if(length(N01) == 1){
        p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
        px <- p.dat[,which(names(p.dat) == p.FE)]
        py2 <- p.dat[,which(names(p.dat) == p.B)]
        # F or E 0.5
        segments(x0 = -1, x1 = N05, y0 = py2[which(px == N05)], y1 = py2[which(px == N05)],
                 col= 'red',lty = 2, lwd=1.5)
        segments(x0 = N05, x1 = N05, y0 = -1, y1 = py2[which(px == N05)],
                 col= 'red',lty = 2, lwd=1.5)
      }
      par(op)
    }

    #Isopleths
    if(type == 'Isopleth'){

      tc_Lc_start <- which(tc_Lc == max(tc_Lc,na.rm=T))
      p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
      px <- p.dat[,which(names(p.dat) == p.FE)]

      if(length(N01) > 1){
        Lc_change <- pes$Lc
        list_Lc_runs <- pes$list_Lc_runs

        Y_R.vec <- list()
        for(fx in 1:length(list_Lc_runs)){
          yr <- list_Lc_runs[[fx]]$Y_R.rel
          Y_R.vec[[fx]] <- yr
        }
        mat_FM_Lc_com.Y <- as.matrix(do.call(cbind,Y_R.vec))
        rownames(mat_FM_Lc_com.Y) <- px
        colnames(mat_FM_Lc_com.Y) <- round(Lc_change/pes$Linf,digits = 2)

        # colours for plot
        pal <- colorRampPalette(rev(c(
          rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1)
        )))

        #plot
        image(x = px,
              y = Lc_change/pes$Linf,
              z = mat_FM_Lc_com.Y, col=pal(100),
              xlab = xlabel1, ylab = 'Lc / Linf')
        contour(x = px,
                y = Lc_change/pes$Linf,
                z = mat_FM_Lc_com.Y, add=TRUE)
        #mtext("Yield", line=0.5, side=3)
      }
    }
  }
}


