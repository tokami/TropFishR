#' @title Plotting catch curve
#
#' @description  This function plots the results from the \code{\link{catchCurve}} model.
#'
#' @param x A list of the class \code{"catchCurve"} containing the results of the
#'      catchCurve model.
#' @param xaxis Character defining if x axis should represent length or age (default: 'age')
#' @param plot_selec logical; if TRUE the regression line is plotted for not fully
#'      exploited length groups and the probability of capture is plotted. This
#'      only works if the \link{catchCurve} was applied with
#'      \code{calc_ogive} == TRUE.
#' @param col a specification for colour of regression points, line and annotation
#' @param cex a numerical value giving the amount by which plotting text and
#'      symbols should be magnified relative to the default.
#' @param xlim limits of x axis
#' @param ylim limits of y axis
#' @param xlab label of x axis. Default display by setting to "default".
#' @param ylab label of y axis. Default display by setting to "default".
#' @param ... standard parameters of plot function
#'
#' @examples
#' \donttest{
#' data(synLFQ3)
#' output <- catchCurve(synLFQ3, calc_ogive = TRUE)
#' plot(output, plot_selec = TRUE)
#' }
#'
#' @details A function to plot the results of the catchCurve model.
#'
#' @importFrom graphics mtext par plot points segments text title
#' @importFrom grDevices dev.cur
#' @importFrom stats fitted
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @method plot catchCurve
#' @export

plot.catchCurve <- function(x, xaxis = 'age', plot_selec = FALSE,
                            col=c('blue',"darkgreen","orange","darkred"),
                            cex = 1.5, xlim = NULL, ylim = NULL, xaxt=NULL,
                            xlab = "default", ylab = "default", ...){
    pes <- x

    ## growth parameters
    if("midLengths" %in% names(pes) | xaxis == "length"){
        if("par" %in% names(pes)){
            Linf <- pes$par$Linf
            K <- pes$par$K
            t0 <- ifelse("t0" %in% names(pes$par), pes$par$t0, 0)
            C <- ifelse("C" %in% names(pes$par), pes$par$C, 0)
            ts <- ifelse("ts" %in% names(pes$par), pes$par$ts, 0)
        }else{
            Linf <- pes$Linf
            K <- pes$K
            t0 <- ifelse("t0" %in% names(pes), pes$t0, 0)
            C <- ifelse("C" %in% names(pes), pes$C, 0)
            ts <- ifelse("ts" %in% names(pes), pes$ts, 0)
        }

        if((is.null(Linf) | is.null(K))) stop(noquote(
                                         "You need to assign values to Linf and K for the catch curve based on length-frequency data!"))
        }


  if(xaxis == 'age'){
    xlabel <- "Age [yrs]"
    if("t_midL" %in% names(pes)){
      xplot <- pes$t_midL
      xlabel <- "Relative age [years - t0]"
      xplotAGE <- pes$t_midL
    }else if("tplusdt_2" %in% names(pes)){
      xplot <- pes$tplusdt_2
    }else if("ln_Linf_L" %in% names(pes)){
      xplot <- pes$ln_Linf_L
      xlabel <- "ln(Linf - L)"
    }else if("classes.num" %in% names(pes) & "midLengths" %in% names(pes)){
        xplot <- pes$classes.num
        xlabel <- "Cohort"
    }else if("classes.num" %in% names(pes)){
        xplot <- pes$classes.num
    }
  }
  if(xaxis == 'length'){
    xplot <- pes$midLengths
    xlabel <- "Length [cm]"
    if("t_midL" %in% names(pes)){
        xplotAGE <- pes$t_midL
    }
  }

  if("lnC_dt" %in% names(pes)){
    yplot <- pes$lnC_dt
    ylabel <- "ln(C / dt)"
    }else if("lnC" %in% names(pes)){
      yplot <- pes$lnC
      ylabel <- "ln(C)"
    }else if("ln_C" %in% names(pes)){
        yplot <- pes$ln_C
        ylabel <- "ln C(L,Linf)"
    }
  if("ln_C" %in% names(pes) & "tplusdt_2" %in% names(pes)) ylabel <- "ln C(t, inf)"


  lm1List <- pes$linear_mod
  Z_lm1List <- pes$Z
  SE_Z_lm1List <- pes$se
  reg_intList <- pes$reg_int
  ## Assumption that Z of smallest selected individuals is most appropriate
  mini <- min(unlist(reg_intList))
  temp <- lapply(reg_intList, function(x) grep(mini,x))
  ind <- sapply(temp, function(x) length(x) > 0)
  cutter <- unlist(reg_intList[ind])


  #for final plot
  minyplot <- ifelse(min(yplot,na.rm=TRUE) < 0, min(yplot,na.rm=TRUE),0)
  maxyplot <- max(yplot,na.rm=TRUE) + 1

  if(is.null(xlim)){
    xlims <- c(min(xplot[which(yplot > 0)], na.rm = TRUE)-0.5,
             max(xplot[which(yplot > 0)], na.rm = TRUE)+0.5)
  }else xlims <- xlim

  if(class(Z_lm1List) == "list"){
      reg_num <- length(Z_lm1List)
  }else{
      reg_num <- 1
  }

  #if(plot_selec & any(names(pes) != "Sest")) writeLines("Please run the catchCurve with calc_ogive == TRUE \nin order to show selectivity plot!")
  if(plot_selec & any(names(pes) == "Sest")){
    maxyplot <- ceiling(pes$intercept)

    if(is.null(ylim)){
      ylims <- c(minyplot, maxyplot)
    }else ylims <- ylim


    if (dev.cur()==1){ # If plot is not open
        opar <- par(mfrow=c(2,1), xpd = FALSE,
                    mar = c(1.2, 4, 1, 1) + 0.1,
                    oma = c(5, 0.5, 1, 2) + 0.1)
        on.exit(par(opar))
    }
    if (dev.cur()==2){ # If plot is open, check if it is a 1x1 plot
        if (all(par()$mfrow == c(2, 1))){
            opar <- par(mfrow=c(2,1), xpd = FALSE,
                        mar = c(1.2, 4, 1, 1) + 0.1,
                        oma = c(5, 0.5, 1, 2) + 0.1)
            on.exit(par(opar))
        }
    }



    ## use user defined labels if given
    if(xlab != "default") xlabel = xlab
    if(ylab != "default") ylabel = ylab


    ## final plot
    plot(x = xplot, y = yplot, ylim = ylims, xaxt=xaxt,
         xlab = '', xaxt = 'n', ylab = ylabel, xlim = xlims,
         cex = cex)

    for(I in 1:reg_num){
          if(reg_num > 1){
              lm1 <- lm1List[[I]]
              reg_int <- reg_intList[[I]]
              Z_lm1 <- Z_lm1List[[I]]
              SE_Z_lm1 <- SE_Z_lm1List[[I]]
          }else{
              if(class(lm1List)=="list"){
                  lm1 <- lm1List[[I]]
              }else{
                  lm1 <- lm1List
              }
              if(class(reg_intList)=="list"){
                  reg_int <- reg_intList[[I]]
              }else{
                  reg_int <- reg_intList
              }
              if(class(Z_lm1List)=="list"){
                  Z_lm1 <- Z_lm1List[[I]]
              }else{
                  Z_lm1 <- Z_lm1List
              }
              if(class(SE_Z_lm1List)=="list"){
                  SE_Z_lm1 <- SE_Z_lm1List[[I]]
              }else{
                  SE_Z_lm1 <- SE_Z_lm1List
              }
          }

          points(x = xplot[reg_int[1]:reg_int[2]], y = yplot[reg_int[1]:reg_int[2]],
                 pch = 19, col = col[I], cex = cex)
          lines(xplot[reg_int[1]:reg_int[2]],fitted(lm1), col=col[I], lwd=1.7)

          if(I == which(ind)){
              temp0 <- seq(0,xplotAGE[reg_int[1]],0.01)
              temp <- predict(lm1,newdata=data.frame(xvar=temp0))
              if(xaxis == 'length'){
                  if("C" %in% names(pes)){
                      temp0 <- VBGF(param = list(Linf = Linf,
                                                 K = K, t0 = t0,
                                                 C = C, ts = ts), t = temp0)
                  }else{
                      ## t0 <- ifelse("t0" %in% names(pes), pes$t0, 0)
                      temp0 <- VBGF(param = list(Linf = Linf, K = K, t0 = t0), t = temp0)
                  }

              }
              lines(temp0,temp, col=col[I], lwd=1.7, lty=2)
          }

          pusr <- par("usr")
          text(x = pusr[2]*0.83, y = pusr[4]-(pusr[4]/(5*I)), labels = paste("Z =",round(Z_lm1,2),"\u00B1",
                                                                round(SE_Z_lm1,2)), col = col[I])


##          mtext(side = 3, line = (reg_num-I+0.3),  text = paste("Z =",round(Z_lm1,2),"+/-",
##                                                                round(SE_Z_lm1,2)), col = col[I])
    }



    plot(pes$Sest ~ xplot, type ='o', xlab = xlabel, xlim = xlims,
         ylab = "Probability of capture")

    if(xaxis == 'length'){
        points(y = 0.5, x=pes$L50,col='red',pch=16)
        segments(x0 = pes$L50, y0 = 0.5, x1 = pes$L50, y1 = 0, col='red',lty=2)
        par(xpd=TRUE)
        title(xlab = xlabel, outer = TRUE, line = 2)
        #text(y=-0.12, x=pes$t50, labels = "t50", col = 'red', xpd=TRUE)
        mtext(text = "L50",side = 1, at = pes$L50,line = 0.3, col = 'red', xpd=TRUE)
    }else{
        points(y = 0.5, x=pes$t50,col='red',pch=16)
        segments(x0 = pes$t50, y0 = 0.5, x1 = pes$t50, y1 = 0, col='red',lty=2)
        par(xpd=TRUE)
        title(xlab = xlabel, outer = TRUE, line = 2)
        #text(y=-0.12, x=pes$t50, labels = "t50", col = 'red', xpd=TRUE)
        mtext(text = "t50",side = 1, at = pes$t50,line = 0.3, col = 'red', xpd=TRUE)
    }




  }else {
      if(is.null(ylim)){
          ylims <- c(minyplot, maxyplot)
      }else ylims <- ylim

      if (dev.cur()==1){ # If plot is not open
          opar <- par(mfrow = c(1,1),
                      mar = c(7, 5, 4, 5) + 0.3)
          on.exit(par(opar))
      }
      if (dev.cur()==2){ # If plot is open, check if it is a 1x1 plot
        if (all(par()$mfrow == c(1, 1))){
            opar <- par(mfrow = c(1,1),
                        mar = c(7, 5, 4, 5) + 0.3)
            on.exit(par(opar))
        }
      }


      ## use user defined labels if given
      if(xlab != "default") xlabel = xlab
      if(ylab != "default") ylabel = ylab

      #final plot
      plot(x = xplot, y = yplot, ylim = ylims,
           xlab = xlabel, ylab = ylabel, xlim = xlims,
           cex = cex)
      par(new=T)

      for(I in 1:reg_num){
          if(reg_num > 1){
              lm1 <- lm1List[[I]]
              reg_int <- reg_intList[[I]]
              Z_lm1 <- Z_lm1List[[I]]
              SE_Z_lm1 <- SE_Z_lm1List[[I]]
          }else{
              if(class(lm1List)=="list"){
                  lm1 <- lm1List[[I]]
              }else{
                  lm1 <- lm1List
              }
              if(class(reg_intList)=="list"){
                  reg_int <- reg_intList[[I]]
              }else{
                  reg_int <- reg_intList
              }
              if(class(Z_lm1List)=="list"){
                  Z_lm1 <- Z_lm1List[[I]]
              }else{
                  Z_lm1 <- Z_lm1List
              }
              if(class(SE_Z_lm1List)=="list"){
                  SE_Z_lm1 <- SE_Z_lm1List[[I]]
              }else{
                  SE_Z_lm1 <- SE_Z_lm1List
              }
          }

          points(x = xplot[reg_int[1]:reg_int[2]], y = yplot[reg_int[1]:reg_int[2]],
                 pch = 19, col = col[I], cex = cex)
          lines(xplot[reg_int[1]:reg_int[2]],fitted(lm1), col=col[I], lwd=1.7)

          pusr <- par("usr")
          text(x = pusr[2]*0.83, y = pusr[4]-(pusr[4]/(12*(1/I))), labels = paste("Z =",round(Z_lm1,2),"\u00B1",
                                                                round(SE_Z_lm1,2)), col = col[I])
          ## mtext(side = 3,line=(reg_num-I+0.3), text = paste("Z =",round(Z_lm1,2),"+/-",
          ##                                                  round(SE_Z_lm1,2)), col = col[I])
    }
  }
}

