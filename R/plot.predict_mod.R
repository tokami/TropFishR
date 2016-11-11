#' @title Plotting prediction models yield per recruit and Thompson & Bell
#'
#' @description This function plots objects of the class "predict_mod", which are results
#'    of the function \code{\link{predict_mod}}.
#'
#' @param x a object of the class 'predict_mod'
#' @param type a character indicating, which type of plot should be displayed in case of
#'    Beverton and Holt's yield per recurit model. Options are either "ypr" (default) for
#'    line plot or "Isopleth" for isopleth plot.
#' @param xaxis1  which x-axis should be plotted? Either "FM" (fishing mortality; default) or "E"
#'    (exploitation rate).
#' @param yaxis1  which (first) y-axis should be plotted? "Y_R" (yield per recruit; default) or
#'    "Y_R.rel" (relative yield per recruit) for type = "ypr". For "Isopleth" in addition: "B_R"
#'    (biomass per recruit) and "B_R.rel" (relative yield per recruit). For Thompson and Bell model
#'    in addition also "value" or "catch" possible.
#' @param yaxis2 which second y-axis should be plotted for type = "ypr"? Either "B_R" (biomass
#'    per recruit; default), "B_R.rel" (relative biomass per recruit), or "B_R.percent"
#'    (percentage biomass per recruit)
#' @param yaxis_iso determines label and scale of y axis of Isopleth graph. Either "Lc"
#'    (default) for length at first capture or "Lc/Linf" for the relation of length
#'    at first capture to the infinite length
#' @param identify logical; indicating whether points in the graph are supposed to be identified by
#'    clicking on them (uses \code{\link{locator}} function). To stop press right mouse click.
#'    (default: TRUE).
#' @param mark logical; if value of choosen points should be displayed in graph (default: TRUE)
#' @param ... optional parameters of plot function
#'
#' @examples
#' \dontrun{
#' # Nemipterus marginatus - age structured data
#' threadfin <- list(Winf = 286,K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4)
#'
#' output <- predict_mod(threadfin, FM_change = seq(0,6,0.1),
#'    tc_change = seq(0.2,1,0.2), type = 'ypr')  #where it is maximal  = MSY
#' plot(output)
#'
#' hake - length structured data
#' data(hake)
#' hake$Lr <- 35
#' select.list <- list(selecType = 'trawl_ogive', L50 = 20, L75 = 24)
#' output <- predict_mod(param = hake, E_change = seq(0,1,0.05),
#'                       curr.E = 0.4, curr.Lc = 40,
#'                       Lc_change = seq(5,80,1), s_list = select.list,
#'                       type = 'ypr', plot = FALSE)
#' plot(output, type = "Isopleth", xaxis1 = "E", yaxis1 = "Y_R.rel",
#'    identify = FALSE)
#'
#'}
#' @importFrom grDevices colorRampPalette dev.new rgb
#' @importFrom graphics mtext par plot axis contour identify image legend lines locator points rect segments text
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.predict_mod <- function(x, type = 'ypr', xaxis1 = "FM",
                             yaxis1 = "Y_R.rel", yaxis2 = "B_R.rel",
                             yaxis_iso = "Lc",
                             identify = FALSE, mark = FALSE, ...){
  pes <- x

  # function for identifying Lc and yield/biomass values in plot
  image.identifier <- function(xyz, markII=TRUE, digits=2){
    intiX <- (xyz$x[2]-xyz$x[1])/2
    intiY <- (xyz$y[2]-xyz$y[1])/2
    newX <- c(xyz$x - intiX,xyz$x[length(xyz$x)]+intiX)
    newY <- c(xyz$y - intiY,xyz$y[length(xyz$y)]+intiY)
    nx <- length(xyz$x)
    ny <- length(xyz$y)
    mesi <- data.frame()
    xy <- locator(1)
    while(!is.null(xy)){
      xbin <- as.numeric(cut(xy$x,newX))
      ybin <- as.numeric(cut(xy$y,newY))
      cat("[",xbin,",",ybin,"] = ",xyz$z[xbin,ybin],"\n",sep='')
      lcc <- xy$y * pes$Linf
      mesi <- rbind(mesi,data.frame(i=xbin,j=ybin,x=xy$x,y=xy$y,z=xyz$z[xbin,ybin],Lc=lcc))
      rm(xy)
      xy <- locator(1)
    }
    if(markII){
      points(mesi$x,mesi$y,pch=19,cex=.5,col="blue")
      text(mesi$x,mesi$y,format(mesi$z,digits=digits),adj=-.2,col="blue")
    }
    colnames(mesi) <- c("i","j",p.FE,"Lc/Linf",p.yield,"Lc")
    # if(mark){
    #   points(xy$x,xy$y,pch=19,cex=.5,col="blue")
    #   text(xy$x,xy$y,format(xyz$z[xbin,ybin],digits=digits),adj=-.2,col="blue")
    # }
    mesi
  }

  #THOMPBELL
  if("totals" %in% names(pes)){
    #save x axis positions
    max_val <- round(max(pes$totV,na.rm=TRUE),digits=0)
    dim_val <- 10 ^ (nchar(max_val)-1)
    max_yiel <- round(max(pes$totY,na.rm=TRUE),digits=0)
    dim_yiel <- 10 ^ (nchar(max_yiel)-1)
    max_bio <- round(max(pes$meanB,na.rm=TRUE),digits=0)
    dim_bio <- 10 ^ (nchar(max_bio)-1)

    #   max_val <- round(max(pes$totV,na.rm=TRUE),digits=0)
    #   dim_val <- 10 ^ (nchar(max_val)-1)

    #op <- par(oma = c(1, 1, 1.5, 1),new=FALSE,mar = c(5, 4, 4, 6) + 0.3)
    plot(pes$Xfact,pes$totV, type ='o',ylab='Value',xlab='F-factor X',
         col ='darkorange', ylim = c(0,ceiling(max_val/dim_val)*dim_val),
         lwd=1.6)
    par(new=TRUE)
    plot(pes$Xfact,pes$totY,type ='o',ylab='',xlab='',
         col='dodgerblue',lwd=1.6,axes=FALSE,
         ylim = c(0,ceiling(max_yiel/dim_yiel) * dim_yiel))
    axis(4,at=pretty(c(0,pes$totY)))
    mtext("Yield", side=4, line=2)
    par(new=TRUE)
    plot(pes$Xfact,pes$meanB,type='o',axes=FALSE,ylab='',xlab='',
         col = 'darkgreen',lwd=1.6,
         ylim = c(0,ceiling(max_bio/dim_bio)*dim_bio))    # draw lines with small intervals: seq(0,max(),0.05) but y as to be dependent of x (formula of calculaiton of y)
    axis(4,at=pretty(c(0,pes$meanB)),line = 3)
    mtext("Biomass", side=4, line=5)

    #par(oma = c(0, 0, 0, 0), new = TRUE)
    legend("top", c("value", "yield", "biomass"), xpd = TRUE,
           horiz = TRUE, inset = c(0, -0.1), bty = "n",lty = 1,seg.len = 0.7,
           col = c('darkorange','dodgerblue','darkgreen'), cex = 0.8,lwd=2,
           text.width=0.3,x.intersp=0.3)
    #par(op)

  }

  # THOMP BELL WITH LC change - ISOPLETHS
  if("mat_FM_Lc_com.C" %in% names(pes)){
    p.yield <- yaxis1
    p.FE <- xaxis1
    p.B <- yaxis2
    xlabel1 <- ifelse(xaxis1 == "FM", "Fishing mortality", "Exploitation rate")
    #ylabel1 <- ifelse(yaxis1 == "Y_R", "Y/R", "rel. Y/R")
    ylabel_iso <- ifelse(yaxis_iso == "Lc", "Lc", "Lc / Linf")

    Lc_change <- pes$Lc_change
    FM_change <- pes$FM_change
    if(p.FE == "FM"){
      px <- FM_change
      if("currents" %in% names(pes)) curr_markX <- pes$currents$curr.F
    }else{
      px <- FM_change/(FM_change + pes$M)
      if("currents" %in% names(pes)) curr_markX <- pes$currents$curr.E
    }

    if(p.yield == "Y_R.rel" | p.yield == "Y_R"){
      pz <- pes$mat_FM_Lc_com.Y
    }else if(p.yield == "B_R.rel" | p.yield == "B_R"){
      pz <- pes$mat_FM_Lc_com.B
    }else if(p.yield == "value"){
      pz <- pes$mat_FM_Lc_com.V
    }else if(p.yield == "catch"){
      pz <- pes$mat_FM_Lc_com.C
    }

    # colours for plot
    pal <- colorRampPalette(rev(c(
      rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1))))

    m <- list(x = px,
              y = Lc_change,
              z = pz)

    if("currents" %in% names(pes)) curr_markY <- pes$currents$curr.Lc

    if(yaxis_iso != "Lc" & "Linf" %in% names(pes)){
      m$y <- Lc_change/pes$Linf
      if("currents" %in% names(pes)) curr_markY <- pes$currents$curr.Lc/pes$Linf
    }

    #plot
    if(identify == TRUE){
      dev.new(noRStudioGD = TRUE)
    }
    image(m, col=pal(100),
          xlab = xlabel1, ylab = ylabel_iso)
    contour(m, add=TRUE)

    if("currents" %in% names(pes) & mark){
      abline(v=curr_markX, col="grey30",lty=2)
      abline(h=curr_markY, col="grey30",lty=2)
    }
    #mtext("Yield", line=0.5, side=3)
    if(identify == TRUE) image.identifier(m)
  }



  # YPR
  if("list_Lc_runs" %in% names(pes) |
     "list_tc_runs" %in% names(pes)){

    # necessary calculations
    FM <- pes$FM
    if("tc" %in% names(pes)) if(!is.null(pes$tc)) tc_Lc <- pes$tc
    if("Lc" %in% names(pes)) if(!is.null(pes$Lc)) tc_Lc <- pes$Lc
    if(is.null(pes$tc) & is.null(pes$Lc)) tc_Lc <- names(pes$list_Lc_runs)
    if("list_tc_runs" %in% names(pes)) list_tc_Lc_runs <- pes$list_tc_runs
    if("list_Lc_runs" %in% names(pes)) list_tc_Lc_runs <- pes$list_Lc_runs

    p.yield <- yaxis1
    p.FE <- xaxis1
    p.B <- yaxis2
    xlabel1 <- ifelse(xaxis1 == "FM", "Fishing mortality", "Exploitation rate")
    ylabel1 <- ifelse(yaxis1 == "Y_R", "Y/R", "rel. Y/R")
    ylabel2 <- ifelse(yaxis2 == "B_R", "B/R", "B/R [%]")
    ylabel_iso <- ifelse(yaxis_iso == "Lc", "Lc", "Lc / Linf")


    df_Es <- pes$df_Es
    if(xaxis1 == "FM"){
      N01 <- df_Es$F01
      N05 <- df_Es$F05
      Nmax <- df_Es$Fmsy
      if(length(N05) == 1){
        legend.lab <- c("F0.1","F0.5","Fmsy")
      }else legend.lab <- c("F0.1","Fmsy")
      if("currents" %in% names(pes)) curr_markX <- pes$currents$curr.F
    }else{
      N01 <- df_Es$E01
      N05 <- df_Es$E05
      Nmax <- df_Es$Emsy
      if(length(N05) == 1){
        legend.lab <- c("E0.1","E0.5","Emsy")
      }else legend.lab <- c("E0.1","Emsy")
      if("currents" %in% names(pes)) curr_markX <- pes$currents$curr.E
    }

    #standard plot (ypr vs. E or FM)
    if(type == 'ypr'){
      #Plot
      #op <- par(mfrow=c(1,1),new=F, mar = c(5, 4, 4, 4) + 0.3)

      #plot Y/R & B/R
      label <- ifelse("tc" %in% names(pes), "tc", "Lc")
      tc_Lc_start <- which(tc_Lc == max(tc_Lc,na.rm=T))
      p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
      py <- p.dat[,which(names(p.dat) == p.yield)]
      px <- p.dat[,which(names(p.dat) == p.FE)]
      offset_text <- py[length(py)] * 0.05
      offset_x <- py[length(px)] * 0.1
      runs <- sapply(strsplit(names(list_tc_Lc_runs),split = "_"), "[[",2)

      plot(px, py, type = 'l', ylim =c(0,max(py, na.rm = TRUE)*1.2),
           ylab = ylabel1, xlab = xlabel1, lty=tc_Lc_start,...)
      text(x = px[length(px)] - offset_x,
           y = py[length(py)] + offset_text,
           labels = bquote(.(label)[.(round(as.numeric(as.character(runs[tc_Lc_start])),2))]))
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
             labels = bquote(.(label)[.(round(as.numeric(as.character(runs[j])),2))]))
      } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area

      # reference points
      if(length(N01) == 1){
        p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
        py <- p.dat[,which(names(p.dat) == p.yield)]
        px <- p.dat[,which(names(p.dat) == p.FE)]
        # F or  E 0.1
        if(N01 == Nmax){
          add_shift <- 1.005
        }else  add_shift <- 1
        segments(x0 = -1, x1 = N01, y0 = py[which(px == N01)],
                 y1 = py[which(px == N01)],
                 col= 'darkgreen',lty = 1, lwd=1.5)
        segments(x0 = N01, x1 = N01, y0 = -1, y1 = py[which(px == N01)],
                 col= 'darkgreen',lty = 1, lwd=1.5)
        # F or E max
        segments(x0 = -1, x1 = Nmax, y0 = py[which(px == Nmax)]*add_shift,
                 y1 = py[which(px == Nmax)]*add_shift,
                 col= 'goldenrod1',lty = 2, lwd=1.5)
        segments(x0 = Nmax*add_shift, x1 = Nmax*add_shift, y0 = -1, y1 = py[which(px == Nmax)],
                 col= 'goldenrod1',lty = 2, lwd=1.5)
        # Legend
        if(length(N05) == 1){
          legend("top", legend = legend.lab, xpd = TRUE, horiz = TRUE,
                 inset = c(0,0), bty = "n", lty = c(1,3,2), col = c("darkgreen","red","goldenrod1"),
                 seg.len = 1,pt.cex = 2, x.intersp = c(0.7,0.7,0.7),merge=TRUE,
                 y.intersp = -2, box.lty=0,cex=0.8, lwd =2)
        }else{
          legend("top", legend = legend.lab, xpd = TRUE, horiz = TRUE,
                 inset = c(0,0), bty = "n", lty = c(1,2), col = c("darkgreen","goldenrod1"),
                 seg.len = 1,pt.cex = 2, x.intersp = c(0.7,0.7,0.7),merge=TRUE,
                 y.intersp = -2, box.lty=0,cex=0.8, lwd =2)
        }


        # current Exploitation rate or fishing mortality
        if(!is.null(pes$currents)){
          currents <- pes$currents
          if(!is.na(currents$curr.E)){
            px <- ifelse(p.FE == "FM",currents$curr.F, currents$curr.E)
            py <- ifelse(p.yield == "Y_R",currents$curr.YPR,currents$curr.YPR.rel)
            points(px,py, pch = 16)
          }
        }
      }

      #Biomass per recruit
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
      if(length(N05) == 1){
        p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
        px <- p.dat[,which(names(p.dat) == p.FE)]
        py2 <- p.dat[,which(names(p.dat) == p.B)]
        # F or E 0.5
        segments(x0 = -1, x1 = N05, y0 = py2[which(px == N05)], y1 = py2[which(px == N05)],
                 col= 'red',lty = 3, lwd=1.5)
        segments(x0 = N05, x1 = N05, y0 = -1, y1 = py2[which(px == N05)],
                 col= 'red',lty = 3, lwd=1.5)
      }
      #par(op)
    }

    #Isopleths
    if(type == 'Isopleth'){

      tc_Lc_start <- which(tc_Lc == max(tc_Lc,na.rm=T))
      p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
      px <- p.dat[,which(names(p.dat) == p.FE)]

      if(length(N01) > 1){
        Lc_change <- pes$Lc
        list_Lc_runs <- pes$list_Lc_runs

        Y_R.vec <- vector("list",length(list_Lc_runs))
        for(fx in 1:length(list_Lc_runs)){
          yr <- list_Lc_runs[[fx]][,which(names(list_Lc_runs[[fx]]) == p.yield)]
          Y_R.vec[[fx]] <- yr
        }
        mat_FM_Lc_com.Y <- as.matrix(do.call(cbind,Y_R.vec))
        rownames(mat_FM_Lc_com.Y) <- round(px,digits=2)
        colnames(mat_FM_Lc_com.Y) <- round(Lc_change/pes$Linf,digits = 2)

        m <- list(x = px,
                  y = Lc_change,
                  z = mat_FM_Lc_com.Y)

        if("currents" %in% names(pes)) curr_markY <- pes$currents$curr.Lc

        if(yaxis_iso != "Lc" & "Linf" %in% names(pes)){
          m$y <- Lc_change/pes$Linf
          if("currents" %in% names(pes)) curr_markY <- pes$currents$curr.Lc/pes$Linf
        }

        # colours for plot
        pal <- colorRampPalette(rev(c(
          rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1)
        )))

        #plot
        if(identify == TRUE){
          dev.new(noRStudioGD = TRUE)
        }
        image(m, col=pal(100),
              xlab = xlabel1, ylab = ylabel_iso)
        contour(m, add=TRUE)
        #mtext("Yield", line=0.5, side=3)
        if("currents" %in% names(pes) & mark){
          abline(v=curr_markX, col="grey30",lty=2)
          abline(h=pes$currents$curr.Lc/pes$Linf, col="grey30",lty=2)
        }
        if(identify == TRUE) image.identifier(m)
      }
    }
  }
}
