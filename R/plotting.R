#' @title Bhattacharya plot
#'
#' @description  This function plots the seperated frequency distributions
#'      and selected regression lines of \link[TropFishR]{Bhattacharya} method.
#'
#' @param x a list of the class \code{"Bhattacharya"} containing the results of
#'      Bhattacharya's method.
#' @param ... additional options of the \code{\link{plot}} function
#' @param analysisPlot logical; indicating wheter the anaylsis graph with the regression
#'      lines should be created
#'
#' @examples
#' \donttest{
#'  data(synLFQ1)
#'  output <- Bhattacharya(param = synLFQ1)
#'  plot(output)
#' }
#'
#' @details This function plots the results of the Bhattacharya method.
#'
#' @importFrom graphics abline hist layout lines par plot points title
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @method plot Bhattacharya
#' @export

plot.Bhattacharya <- function(x, analysisPlot = TRUE, ...){
    pes <- x

    s.l.df <- pes$Lmean_SD_list
    a.b.df <- pes$regression_lines
    bhat.results <- pes$bhat_results
    temp.list <- pes$distributions


    midLengths <- bhat.results$mean.length.classes
    N <- bhat.results$N1.plus
    frequis <- rep(midLengths, N)

    colour.xy <- c('blue','darkgreen','red','goldenrod2','purple','orange',
                   'lightgreen','skyblue','brown','darkblue','darkorange','darkred')

    h <- hist(frequis, breaks = 50, plot = FALSE)

    ## histogram
    maxTemps <- rep(NA,length(temp.list))
    for(maxis in 1:length(temp.list)){
        maxTemps[maxis] <- max(temp.list[[maxis]]$yyfit, na.rm = TRUE)
    }

    maxyyfit <- max(maxTemps, na.rm = TRUE)
    maxlim <- ifelse(maxyyfit > max(h$density), maxyyfit, max(h$density))

    if(analysisPlot == FALSE){
        par(mfrow = c(1,1))

        hist(frequis, breaks = 50, main = "", xlab = "L",
             probability = TRUE, ylim = c(0,maxlim))

        for(cohorti in 1:length(temp.list)){
            lines(temp.list[[cohorti]]$xfit, temp.list[[cohorti]]$yyfit,
                  col=colour.xy[cohorti], lwd=1.5)
        }
    }

    if(analysisPlot == TRUE){
        par(mfrow = c(2,1),
            oma = c(6.5, 0.5, 2, 1) + 0.1,
            mar = c(0, 4, 0.5, 0) + 0.1)
        layout(matrix(c(1,2),nrow=2), heights = c(1,2.5))

        hist(frequis, breaks = 50, main = "", xaxt = 'n',
             probability = TRUE, ylim = c(0,maxlim))

        for(cohorti in 1:length(temp.list)){
            lines(temp.list[[cohorti]]$xfit, temp.list[[cohorti]]$yyfit,
                  col=colour.xy[cohorti], lwd=1.5)
        }

        plot(bhat.results$L, bhat.results$delta.log.N1.plus, pch = 4,
             col = 'black', xlab = "",...,
             ylab = expression(paste(delta," log N+")))
        abline(h = 0)
        seqi <- seq(4,length(names(bhat.results)),by = 7) ## depending on order of bhat.results table

        for(coli in 1:length(seqi)){
            sta <- a.b.df[coli,4]
            end <- a.b.df[coli,5]
            points(bhat.results$L[sta:end],
                   bhat.results[sta:end,seqi[coli]],
                   col = colour.xy[coli], pch = 16)
        }
        for(coli in 1:length(seqi)){
            ai <- a.b.df[coli,2]
            bi <- a.b.df[coli,3]
            abline(a = ai, b=bi, col=colour.xy[coli], lwd = 1.5)
        }
        title(xlab = "L", outer = TRUE, line = 2.5)
    }
}

#' @title Plotting catch curve
#'
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
                            cex = 1.5, xlim = NULL, ylim = NULL,
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
            xlabel <- "Relative age [yrs]"
            xplotAGE <- pes$t_midL
        }else if("tplusdt_2" %in% names(pes)){
            xplot <- pes$tplusdt_2
        }else if("ln_Linf_L" %in% names(pes)){
            xplot <- pes$ln_Linf_L
            xlabel <- "ln(Linf - L)"
        }else if("classes.num" %in% names(pes)) xplot <- pes$classes.num
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
    Z_lm1List <- pes$par$Z
    SE_Z_lm1List <- pes$Z_se
    reg_intList <- pes$reg_int
    ## Assumption that Z of smallest selected individuals is most appropriate
    mini <- min(unlist(reg_intList))
    temp <- lapply(reg_intList, function(x) grep(mini,x))
    ind <- sapply(temp, function(x) length(x) > 0)
    cutter <- unlist(reg_intList[ind])


    ##for final plot
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

    ##if(plot_selec & any(names(pes) != "Sest")) writeLines("Please run the catchCurve with calc_ogive == TRUE \nin order to show selectivity plot!")
    if(plot_selec & any(names(pes) == "Sest")){
        maxyplot <- ceiling(pes$intercept)

        if(is.null(ylim)){
            ylims <- c(minyplot, maxyplot)
        }else ylims <- ylim


        if (dev.cur()==1){ ## If plot is not open
            opar <- par(mfrow=c(2,1), xpd = FALSE,
                        mar = c(1.2, 4, 1, 1) + 0.1,
                        oma = c(5, 0.5, 1, 2) + 0.1)
            on.exit(par(opar))
        }
        if (dev.cur()==2){ ## If plot is open, check if it is a 1x1 plot
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
        plot(x = xplot, y = yplot, ylim = ylims,
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
                        temp0 <- VBGF(pars = list(Linf = Linf,
                                                   K = K, t0 = t0,
                                                   C = C, ts = ts), t = temp0)
                    }else{
                        ## t0 <- ifelse("t0" %in% names(pes), pes$t0, 0)
                        temp0 <- VBGF(pars = list(Linf = Linf, K = K, t0 = t0), t = temp0)
                    }

                }
                lines(temp0,temp, col=col[I], lwd=1.7, lty=2)
            }

            pusr <- par("usr")
            text(x = pusr[2]*0.85, y = pusr[4]-(pusr[4]/(5*I)), labels = paste("Z =",round(Z_lm1,2),"+/-",
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
            ##text(y=-0.12, x=pes$t50, labels = "t50", col = 'red', xpd=TRUE)
            mtext(text = "L50",side = 1, at = pes$L50,line = 0.3, col = 'red', xpd=TRUE)
        }else{
            points(y = 0.5, x=pes$t50,col='red',pch=16)
            segments(x0 = pes$t50, y0 = 0.5, x1 = pes$t50, y1 = 0, col='red',lty=2)
            par(xpd=TRUE)
            title(xlab = xlabel, outer = TRUE, line = 2)
            ##text(y=-0.12, x=pes$t50, labels = "t50", col = 'red', xpd=TRUE)
            mtext(text = "t50",side = 1, at = pes$t50,line = 0.3, col = 'red', xpd=TRUE)
        }




    }else {
        if(is.null(ylim)){
            ylims <- c(minyplot, maxyplot)
        }else ylims <- ylim

        if (dev.cur()==1){ ## If plot is not open
            opar <- par(mfrow = c(1,1),
                        mar = c(7, 5, 4, 5) + 0.3)
            on.exit(par(opar))
        }
        if (dev.cur()==2){ ## If plot is open, check if it is a 1x1 plot
            if (all(par()$mfrow == c(1, 1))){
                opar <- par(mfrow = c(1,1),
                            mar = c(7, 5, 4, 5) + 0.3)
                on.exit(par(opar))
            }
        }


        ## use user defined labels if given
        if(xlab != "default") xlabel = xlab
        if(ylab != "default") ylabel = ylab

        ##final plot
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
            text(x = pusr[2]*0.85, y = pusr[4]-(pusr[4]/(12*(1/I))), labels = paste("Z =",round(Z_lm1,2),"+/-",
                                                                                    round(SE_Z_lm1,2)), col = col[I])
            ## mtext(side = 3,line=(reg_num-I+0.3), text = paste("Z =",round(Z_lm1,2),"+/-",
            ##                                                  round(SE_Z_lm1,2)), col = col[I])
        }
    }
}


#' @title Plotting of length frequency data (with VBGF curves)
#'
#' @description This function plots length frequency (lfq) samples sequentially
#'  arranged in time. An object
#'  of "lfq" class is obatined by applying the \code{\link{lfqRestructure}}
#'  function. In case growth
#'  parameters are known, von Bertalanffy growth curves can be plotted through the lfq samples.
#'
#' @param x a list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column);
#' }
#' @param Fname indicating whether restructured ("rcounts") or original frequencies ("catch") should be
#'    displayed (default: "rcounts")
#' @param par a list with following growth parameters (default NULL):
#'  \itemize{
#'   \item \strong{Linf} asymptotic length,
#'   \item \strong{K} growth coefficient,
#'   \item \strong{ta} time at length zero,
#'   \item \strong{C} amplitude of growth oscillation (optional),
#'   \item \strong{ts} summer point (optional);
#' }
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param rel logical; defines if relative numbers per length class should be plotted (relative to
#'   the sample size per sampling time, e.g. month). Default: FALSE.
#' @param y an optional second list of class "lfq" consisting of same parameters as x. This allows to plot
#'   samples from different sources (e.g. different fleets) on top of eachother, classes and dates have to
#'   correspond at least partially. Default is NA.
#' @param curve.col colour of growth curves (default: 1)
#' @param hist.sc defines the scaling factor to use for maximum histogram extent (x-axis
#'    direction). The default setting of hist.sc=0.5 will result in a maximum distance equal
#'    to half the distance between closest sample dates (i.e. ensures no overlap and full
#'    plotting within the plot region).
#' @param hist.col vector of 2 values defining coloring to use on negative and positive
#'    histogram bars (default: hist.col=c("white", "black", "orange", "darkgreen"))
#' @param image.col colour of image, by default (NULL) red and blue colours
#'    are used. To remove image coloring, set image.col=NA.
#' @param region.col colour of plotting region. Will overwrite image.col
#' (default: region.col=NULL)
#' @param zlim the minimum and maximum z values for which colors should be
#'    plotted (default : NULL).
#' @param zlimtype indicating if zlim should be based on the range of
#'    the catch ("range") or based on the maximum absolute value in Fname
#'    ("balanced", default). This parameter is only considered if zlim is NULL.
#' @param date.axis the style of the x axis. By default the "traditional"
#'    approach is used with years under the months. Alternatively, by using
#'    "modern" the date is plotted in one line according to the chosen
#'     format \code{date.format}.
#' @param date.at the points at which tick-marks are to be drawn. Non-finite
#'    (infinite, NaN or NA) values are omitted. By default it is
#'    seq(as.Date("1500-01-01"), as.Date("2500-01-01"), by="months")
#' @param date.format format of date  (default : "\%y-\%b")
#' @param xlab label of x axis (default : "")
#' @param ylab label of y axis (default : "Length classes")
#' @param draw logical; indicating whether growth curves should be added to
#'     lfq plot if parameters are provided (default : TRUE)
#' @param ... additional options of the plot function
#'
#'
#' @examples
#' data(alba)
#' res <- lfqRestructure(alba)
#'
#' # simple plot or reconstructed frequencies
#' plot(x = res, Fname = "rcounts")
#'
#' # add VBGF curves
#' plot(res, Fname = "rcounts", par = list(Linf = 14, K = 1.1, ta = 0.3))
#'
#' # add soVBGF curves, adjust hist.sc and xlim
#' plot(res, Fname = "catch", curve.col=4,
#'   par = list(Linf = 14, K = 1.1, ta = 0.3, C = 0.2, ts = 0.75),
#'   hist.sc = 0.9,
#'   xlim=range(res$dates)+c(-30, 0)
#' )
#'
#' # adjust image colors
#' plot(res, Fname = "rcounts", image.col = NA )
#' plot(res, Fname = "rcounts", image.col = rev(cm.colors(21)) )
#' plot(res, Fname = "rcounts", image.col = colorRampPalette(c("red","grey90","green"))(21))
#'
#' # solid plot region color
#' plot(res, xlim=range(res$dates)+c(-60, 60),
#'   hist.sc=0.75, image.col="grey90") # leaves gaps
#' plot(res, xlim=range(res$dates)+c(-60, 60),
#'   hist.sc=0.75, region.col="grey90") # full coverage
#'
#' # low-level plot additions
#' plot(res)
#' abline(h=4, lty=2)
#' mtext("Restructured frequencies (MA=5)", line=0.25, side=3)
#'
#'
#'
#' @details This function uses \code{\link{lfqFitCurves}} when growth
#'    parameters are provided to plot growth curves, this can be turned off with
#'    \code{draw} = FALSE.
#'
#' @importFrom grDevices colorRampPalette dev.new dev.off recordPlot rgb
#' @importFrom graphics abline axis box lines mtext par plot rect polygon axis.Date
#' @importFrom stats update
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a revised and expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. 1981. The relationship between gill surface area and growth performance in fish:
#' a generalization of von Bertalanffy's theory of growth. \emph{Meeresforschung}. 28:205-211
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Pauly, D., 1985. On improving operation and use of ELEFAN programs. Part I: Avoiding
#' "drift" of K towards low values. \emph{ICLARM Conf. Proc.}, 13-14
#'
#' Pauly, D., 1987. A review of the ELEFAN system for analysis of length-frequency data in
#' fish and aquatic invertebrates. \emph{ICLARM Conf. Proc.}, (13):7-34
#'
#' Pauly, D. and G. R. Morgan (Eds.), 1987. Length-based methods in fisheries research.
#' (No. 13). WorldFish
#'
#' Pauly, D. and G. Gaschuetz. 1979. A simple method for fitting oscillating length growth data, with a
#' program for pocket calculators. I.C.E.S. CM 1979/6:24. Demersal Fish Comittee, 26 p.
#'
#' Pauly, D. 1984. Fish population dynamics in tropical waters: a manual for use with programmable
#' calculators (Vol. 8). WorldFish.
#'
#' Quenouille, M. H., 1956. Notes on bias in estimation. \emph{Biometrika}, 43:353-360
#'
#' Somers, I. F., 1988. On a seasonally oscillating growth function. ICLARM Fishbyte 6(1): 8-11.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2): 407 p.
#'
#' Tukey, J., 1958. Bias and confidence in not quite large samples.
#' \emph{Annals of Mathematical Statistics}, 29: 614
#'
#' Tukey, J., 1986. The future of processes of data analysis. In L. V. Jones (Eds.),
#' The Collected Works of John W. Tukey philosophy and principles of data analysis:
#' 1965-1986 (Vol. 4, pp. 517-549). Monterey, CA, USA: Wadsworth & Brooks/Cole
#'
#' @method plot lfq
#' @export

plot.lfq <- function(x,
                     Fname = "rcounts",  ## alternative : "catch"
                     par = NULL,
                     agemax = NULL,
                     rel = FALSE,
                     y = NA,
                     curve.col = 1,
                     hist.sc = 0.5,
                     hist.col = c("white", "black", "orange","darkgreen"),
                     image.col = NULL,
                     region.col = NULL,
                     zlim = NULL,
                     zlimtype = "balanced",   ## alternative : "range"
                     date.axis = "traditional",  ## alternative : "modern"
                     date.at = seq(as.Date("1500-01-01"), as.Date("2500-01-01"), by="months"),
                     date.format = "'%y-%b", xlab = "", ylab = "Length classes",
                     draw = TRUE,
                     ...){

    dates <- x$dates
    classes <- x$midLengths
    catch <- get(Fname, x)

    ## combine lfq data sets (e.g. different fleets)
    if(any(!is.na(y))){
        datesY <- y$dates
        classesY <- y$midLengths
        catchY <- get(Fname, y)

        if(diff(classes)[1] != diff(classesY)[1]) stop("The bin sizes do not fit eachother")
        mergi <- merge(data.frame(dates=dates,x=dates),
                       data.frame(dates=datesY,y=datesY),
                       by="dates",all=TRUE)
        mergi2 <- merge(data.frame(classes=classes,x=classes),
                        data.frame(classes=classesY,y=classesY),
                        by="classes",all=TRUE)
        indY <- which(is.na(mergi2$y) & mergi2$classes > max(classesY))
        matY <- matrix(0, nrow=length(indY), ncol=ncol(catchY))
        catchY <- rbind(catchY,matY)
        indY <- which(is.na(mergi2$y) & mergi2$classes < min(classesY))
        matY <- matrix(0, nrow=length(indY), ncol=ncol(catchY))
        catchY <- rbind(matY,catchY)
        ind <- which(is.na(mergi2$x) & mergi2$classes > max(classes))
        mat <- matrix(0, nrow=length(ind), ncol=ncol(catch))
        catch <- rbind(catch,mat)
        ind <- which(is.na(mergi2$x) & mergi2$classes < min(classes))
        mat <- matrix(0, nrow=length(ind), ncol=ncol(catch))
        catch <- rbind(mat,catch)

        designMat <- matrix(0, ncol=length(mergi$dates), nrow=length(mergi2$classes))
        temp <- designMat
        ind = 1
        for(i in which(!is.na(mergi$x))){
            temp[,i] <- catch[,ind]
            ind <- ind + 1
        }
        catch <- temp
        temp <- designMat
        ind = 1
        for(i in which(!is.na(mergi$y))){
            temp[,i] <- catchY[,ind]
            ind <- ind + 1
        }
        catchY <- temp
    }

    ## display relative catches (relative to number of samples per month)
    if(rel){
        catchrel <- catch
        for(i in 1:ncol(catch)){
            catchrel[,i] <- catch[,i]/colSums(catch, na.rm = TRUE)[i]
        }
        catch <- catchrel
        catch[is.nan(catch)] <- 0

        ## combined lfq plot
        if(any(!is.na(y))){
            catchrelY <- catchY
            for(i in 1:ncol(catchY)){
                catchrelY[,i] <- catchY[,i]/colSums(catchY, na.rm = TRUE)[i]
            }
            catchY <- catchrelY
            catchY[is.nan(catchY)] <- 0
        }
    }


    ## bin height scaling
    sc <- unclass(min(diff(dates)) * hist.sc / max(abs(catch)))

    if(any(!is.na(y))){
        ## bin height scaling
        scY <- unclass(min(diff(dates)) * hist.sc / max(abs(catchY)))
    }

    bin.width <- diff(classes)
    bin.lower <- classes - c(bin.width[1], bin.width)/2
    bin.upper <- classes + c(bin.width, bin.width[length(bin.width)])/2

    ## image colour
    if(is.null(image.col)){
        pal <- colorRampPalette(c(rgb(1,0.8,0.8), rgb(1,1,1), rgb(0.8,0.8,1)))
        image.col <- pal(21)
    }
    if(!is.null(region.col)){
        image.col <- NA
    }

    ## zlim value + type
    if(is.null(zlim) & zlimtype == "balanced"){
        zlim = c(-1,1) * max(abs(catch), na.rm=TRUE)
    }
    if(is.null(zlim) & zlimtype == "range"){
        zlim = range(catch, na.rm = TRUE)
    }

    ## Initial plot
    if(any(!is.na(y))){
        catchComb <- catch + catchY
        image(
            x=mergi$dates, y=mergi2$classes, z=t(catchComb), col=image.col, zlim=zlim,
            xaxt="n", xlab = xlab, ylab = ylab, ...)
    }else{
        image(
            x=dates, y=classes, z=t(catch), col=image.col, zlim=zlim,
            xaxt="n", xlab = xlab, ylab = ylab, ...)
    }

    if(!is.null(region.col)){
        usr <- par()$usr
        if(par()$xlog) usr[1:2] <- 10^usr[1:2]
        if(par()$ylog) usr[3:4] <- 10^usr[3:4]
        rect(usr[1], usr[3], usr[2], usr[4], col=region.col)
    }

    ## add time axis
    if(date.axis == "modern"){
        axis.Date(side = 1, x=dates, at=date.at, format = date.format)
    }else if(date.axis == "traditional"){
        axis.Date(side = 1, x = dates, at = date.at, format = "%b")
        year <- seq(min(as.numeric(format(dates, "%Y"))), max(as.numeric(format(dates, "%Y"))), 1)
        date_seq <- seq.Date(dates[1],dates[length(dates)], by = "month")
        date_label <- format(date_seq, "%m")
        year_pre <- which(date_label %in% "01")
        if(!(1 %in% year_pre)) year_pre <- c(1,which(date_label %in% "01"))
        dates_for_years <- as.Date(paste(format(date_seq,"%Y"),date_label,"01",sep="-"))
        year_ticks <- dates_for_years[year_pre]
        mtext(side = 1, at = year_ticks, text = year, line = 2.5)
    }

    ## Histograms
    if(any(!is.na(y))){
        bin.width <- diff(mergi2$classes)
        bin.lower <- mergi2$classes - c(bin.width[1], bin.width)/2
        bin.upper <- mergi2$classes + c(bin.width, bin.width[length(bin.width)])/2
        for(i in seq(length(mergi$dates))){
            score.sc <- catch[,i] * sc
            score.scY <- catchY[,i] * scY
            for(j in seq(mergi2$classes)){
                polygon(
                    x = c(mergi$dates[i], mergi$dates[i], mergi$dates[i]-score.sc[j], mergi$dates[i]-score.sc[j]),
                    y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
                    col = hist.col[(score.sc[j]>0)+1],
                    border = "grey20", lwd = 1)
                polygon(
                    x = c(mergi$dates[i], mergi$dates[i], mergi$dates[i]-score.scY[j], mergi$dates[i]-score.scY[j]),
                    y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
                    col = hist.col[(score.scY[j]>0)+3],
                    border = "grey20", lwd = 1)

            }
        }
    }else{
        for(i in seq(length(dates))){
            score.sc <- catch[,i] * sc
            for(j in seq(classes)){
                ## if(score.sc[j] != 0){
                polygon(
                    x = c(dates[i], dates[i], dates[i]-score.sc[j], dates[i]-score.sc[j]),
                    y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
                    col = hist.col[(score.sc[j]>0)+1],
                    border = "grey20", lwd = 1)
                ## }
            }
        }
    }

    ## optional addition of cohort growth curves
    if("par" %in% names(x) & is.null(par) & draw){
        Lt <- lfqFitCurves(lfq = x, par = x$par,
                           agemax = x$agemax, draw = TRUE, col=curve.col
                           )
    }
    if(!is.null(par) & draw){
        Lt <- lfqFitCurves(x, par = par,
                           agemax = agemax, draw = TRUE, col=curve.col
                           )
    }

    ## frame
    box()
}

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
#' @param contour used in combination with the Isopleth graph. Usage
#'    can be logical (e.g. TRUE) or providing a numeric which indicates the
#'    number of levels (\code{nlevels} in \code{\link{contour}}). By default TRUE.
#' @param xlab Label of x-axis. If set to NA, then default "Fishing mortality", or
#'    "Exploitation rate" is used.
#' @param ylab1 Label of y-axis. If set to NA, then default "Yield" is used for the
#'     Thompson and Bell model, "Lc" is used for the Isopleth graph, and "Y/R" is used for ypr.
#' @param ylab2 Label of second y-axis. If set to NA, then default "Biomass" is used for
#'    the Thompson and Bell model and "B/R" for ypr.
#' @param ylab3 Label of third y-axis. If set to NA, then default "Value" is used
#'    for the Thompson and Bell model.
#' @param ... optional parameters of plot function
#'
#' @examples
#' \dontrun{
#' # Nemipterus marginatus - age structured data
#' threadfin <- list(par = list(Winf = 286,K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4))
#'
#' output <- predict_mod(threadfin, FM_change = seq(0,6,0.1),
#'    tc_change = seq(0.2,1,0.2), type = 'ypr')
#' plot(output)
#'
#' # hake - length structured data
#' data(hake)
#' hake$par$Lr <- 35
#' select.list <- list(selecType = 'trawl_ogive', L50 = 20, L75 = 24)
#' output <- predict_mod(lfq = hake, FM_change = seq(0,0.4,0.05),
#'                       curr.E = 0.4, curr.Lc = 40,
#'                       Lc_change = seq(5,80,1), s_list = select.list,
#'                       type = 'ypr', plot = FALSE)
#' plot(output, type = "Isopleth", xaxis1 = "FM", yaxis1 = "Y_R.rel",
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
#' @method plot predict_mod
#' @export

plot.predict_mod <- function(x, type = 'ypr', xaxis1 = "FM",
                             yaxis1 = "Y_R.rel", yaxis2 = "B_R.rel",
                             yaxis_iso = "Lc",
                             identify = FALSE, mark = FALSE, contour = TRUE,
                             xlab = NA, ylab1 = NA, ylab2 = NA, ylab3 = NA,
                             ...){
    pes <- x

    ## function for identifying Lc and yield/biomass values in plot
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
        ## if(mark){
        ##   points(xy$x,xy$y,pch=19,cex=.5,col="blue")
        ##   text(xy$x,xy$y,format(xyz$z[xbin,ybin],digits=digits),adj=-.2,col="blue")
        ## }
        mesi
    }

    ##THOMPBELL
    if("totY" %in% names(pes)){
        df_Es <- pes$df_Es

        if(xaxis1 == "FM"){
            px <- pes$FM_change
            xlabel1 <- "Fishing mortality"
            if(pes$FM_relative){
                xlabel1 <- "rel. Fishing mortality"
            }

            N05 <- df_Es$F05
            Nmax <- df_Es$Fmax
            if(length(N05) == 1){
                legend.lab <- c("F0.5","Fmax")
            }else legend.lab <- c("Fmax")
            if("currents" %in% names(pes)) curr_markX <- pes$currents$curr.F
        }else{
            px <- pes$E_change
            xlabel1 <- "Exploitation rate"
            if(pes$FM_relative){
                xlabel1 <- "rel. Exploitation rate"
            }

            N05 <- df_Es$E05
            Nmax <- df_Es$Emax
            if(length(N05) == 1){
                legend.lab <- c("E0.5","Emax")
            }else legend.lab <- c("Emax")
            if("currents" %in% names(pes)) curr_markX <- pes$currents$curr.E
        }
        if(!is.na(xlab[1])){
            xlabel1 <- xlab
        }
        if(is.na(ylab1[1])){
            ylabel1 <- "Yield"
        }
        if(!is.na(ylab1[1])){
            ylabel1 <- ylab1
        }
        if(is.na(ylab2[1])){
            ylabel2 <- "Biomass"
        }
        if(!is.na(ylab2[1])){
            ylabel2 <- ylab2
        }
        if(is.na(ylab3[1])){
            ylabel3 <- "Value"
        }
        if(!is.na(ylab3[1])){
            ylabel3 <- ylab3
        }



        ##save x axis positions
        max_val <- round(max(pes$totV,na.rm=TRUE),digits=0)
        dim_val <- 10 ^ (nchar(max_val)-1)
        max_yiel <- round(max(pes$totY,na.rm=TRUE),digits=0)
        dim_yiel <- 10 ^ (nchar(max_yiel)-1)
        max_bio <- round(max(pes$meanB,na.rm=TRUE),digits=0)
        dim_bio <- 10 ^ (nchar(max_bio)-1)

        ##   max_val <- round(max(pes$totV,na.rm=TRUE),digits=0)
        ##   dim_val <- 10 ^ (nchar(max_val)-1)


        py <- pes$totY[1:length(px)]
        py2 <- pes$meanB[1:length(px)]

        ##op <- par(oma = c(1, 1, 1.5, 1),new=FALSE,mar = c(5, 4, 4, 6) + 0.3)
        plot(px,py, type ='l',ylab=ylabel1,xlab= xlabel1,
             col ='black', ylim = c(0,ceiling(max_yiel/dim_yiel)*dim_yiel),
             lwd=1.6)
        ## F or E max
        segments(x0 = -1, x1 = Nmax, y0 = py[which(px == Nmax)],
                 y1 = py[which(px == Nmax)],
                 col= 'goldenrod1',lty = 2, lwd=1.6)
        segments(x0 = Nmax, x1 = Nmax, y0 = -1, y1 = py[which(px == Nmax)],
                 col= 'goldenrod1',lty = 2, lwd=1.6)

        ## current Exploitation rate or fishing mortality
        if(!is.null(pes$currents) & mark){
            currents <- pes$currents
            if(!is.na(currents$curr.E) & yaxis1 == "Y_R" | yaxis1 == "Y_R.rel"){
                px1 <- ifelse(xaxis1 == "FM",currents$curr.F, currents$curr.E)
                py1 <- currents$curr.Y
                points(px1,py1, pch = 16, col="grey30")
                abline(v=px1, col="grey30",lty=2)
            }
        }

        ## Biomass
        par(new=TRUE)
        plot(px,py2,type ='l',ylab='',xlab='',
             col='blue',lwd=1.6,axes=FALSE)
        axis(4,at=pretty(c(0,max(pes$meanB))),col = "blue", col.axis="blue")
        mtext(ylabel2, side=4, line=2.5, col = "blue", cex=1)
        ## F or E 05
        segments(x0 = -1, x1 = N05, y0 = py2[which(px == N05)], y1 = py2[which(px == N05)],
                 col= 'red',lty = 3, lwd=1.5)
        segments(x0 = N05, x1 = N05, y0 = -1, y1 = py2[which(px == N05)],
                 col= 'red',lty = 3, lwd=1.5)
        ## current Exploitation rate or fishing mortality
        if(!is.null(pes$currents) & mark){
            currents <- pes$currents
            if(!is.na(currents$curr.E) & yaxis1 == "B_R" | yaxis1 == "B_R.rel"){
                px1 <- ifelse(xaxis1 == "FM",currents$curr.F, currents$curr.E)
                py1 <- currents$curr.B
                points(px1,py1, pch = 16, col="grey30")
                abline(v=px1, col="grey30",lty=2)
            }
        }
        ## Legend
        legend("top", legend = legend.lab, xpd = TRUE, horiz = TRUE,
               inset = c(0,0), bty = "n", lty = c(1,2), col = c("red","goldenrod1"),
               seg.len = 1,pt.cex = 2, x.intersp = c(0.7,0.7),merge=TRUE,
               y.intersp = -2, box.lty=0,cex=0.8, lwd =2)

        ## Value if present
        if(any(pes$totV != 0)){
            py3 <- pes$totV[1:length(px)]
            par(new=TRUE)
            plot(px,py3,type='l',axes=FALSE,ylab='',xlab='',
                 col = 'darkgreen',lwd=1.6,
                 ylim = c(0,ceiling(max_val/dim_val)*dim_val))    ## draw lines with small intervals: seq(0,max(),0.05) but y as to be dependent of x (formula of calculaiton of y)
            axis(4,at=pretty(c(0,pes$totV)),line = 3.6,col.axis="darkgreen",col="darkgreen")
            mtext(ylabel3, side=4, line=5.7,col="darkgreen")
        }

        ## ##par(oma = c(0, 0, 0, 0), new = TRUE)
        ## legend("top", c("value", "yield", "biomass"), xpd = TRUE,
        ##        horiz = TRUE, inset = c(0, -0.1), bty = "n",lty = 1,seg.len = 0.7,
        ##        col = c('darkorange','dodgerblue','darkgreen'), cex = 0.8,lwd=2,
        ##        text.width=0.3,x.intersp=0.3)
        ## ##par(op)
    }

    ## THOMP BELL WITH LC change - ISOPLETHS
    if("mat_FM_Lc_com.C" %in% names(pes)){
        p.yield <- yaxis1
        p.FE <- xaxis1
        p.B <- yaxis2
        xlabel1 <- ifelse(xaxis1 == "FM", "Fishing mortality", "Exploitation rate")
        if(pes$FM_relative){
            xlabel1 <- ifelse(xaxis1 == "FM", "rel. Fishing mortality", "rel. Exploitation rate")
        }
        ##ylabel1 <- ifelse(yaxis1 == "Y_R", "Y/R", "rel. Y/R")
        ylabel_iso <- ifelse(yaxis_iso == "Lc", "Lc", "Lc / Linf")
        if(!is.na(xlab[1])){
            xlabel1 <- xlab
        }
        if(!is.na(ylab1[1])){
            ylabel_iso <- ylab1
        }

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

        ## colours for plot
        pal <- colorRampPalette(rev(c(
            rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1))))
        if(yaxis1 == "B_R" | yaxis1 == "B_R.rel"){
            pal <- colorRampPalette(c(
                rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1)))
        }

        m <- list(x = px,
                  y = Lc_change,
                  z = pz)

        if("currents" %in% names(pes)) curr_markY <- pes$currents$curr.Lc

        if(yaxis_iso != "Lc" & "Linf" %in% names(pes)){
            m$y <- Lc_change/pes$Linf
            if("currents" %in% names(pes)) curr_markY <- pes$currents$curr.Lc/pes$Linf
        }        

        ##plot
        if(identify == TRUE){
            dev.new(noRStudioGD = TRUE)
        }
        image(m, col=pal(100),
              xlab = xlabel1, ylab = ylabel_iso)
        if(is.numeric(contour)){
            contour(m, add=TRUE, nlevels = contour)
        }else if(contour == TRUE){
            contour(m, add=TRUE)
        }

        if("currents" %in% names(pes) & mark){
            points(x = curr_markX, y=curr_markY, pch=16, col="grey30")
            abline(v=curr_markX, col="grey30",lty=2)
            abline(h=curr_markY, col="grey30",lty=2)
        }
        ##mtext("Yield", line=0.5, side=3)
        if(identify == TRUE) image.identifier(m)
    }



    ## YPR
    if("list_Lc_runs" %in% names(pes) |
       "list_tc_runs" %in% names(pes)){

        ## necessary calculations
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
        if(pes$FM_relative){
            xlabel1 <- ifelse(xaxis1 == "FM", "rel. Fishing mortality", "rel. Exploitation rate")
        }
        ylabel1 <- ifelse(yaxis1 == "Y_R", "Y/R", "rel. Y/R")
        ylabel2 <- ifelse(yaxis2 == "B_R", "B/R", "B/R [%]")
        ylabel_iso <- ifelse(yaxis_iso == "Lc", "Lc", "Lc / Linf")
        if(!is.na(xlab[1])){
            xlabel1 <- xlab
        }
        if(!is.na(ylab1[1])){
            ylabel1 <- ylab1
        }
        if(!is.na(ylab2[1])){
            ylabel2 <- ylab2
        }
        if(!is.na(ylab1[1])){
            ylabel_iso <- ylab1
        }



        df_Es <- pes$df_Es
        if(xaxis1 == "FM"){
            N01 <- df_Es$F01
            N05 <- df_Es$F05
            Nmax <- df_Es$Fmax
            if(length(N05) == 1){
                legend.lab <- c("F0.1","F0.5","Fmax")
            }else legend.lab <- c("F0.1","Fmax")
            if("currents" %in% names(pes)) curr_markX <- pes$currents$curr.F
        }else{
            N01 <- df_Es$E01
            N05 <- df_Es$E05
            Nmax <- df_Es$Emax
            if(length(N05) == 1){
                legend.lab <- c("E0.1","E0.5","Emax")
            }else legend.lab <- c("E0.1","Emax")
            if("currents" %in% names(pes)) curr_markX <- pes$currents$curr.E
        }

        ##standard plot (ypr vs. E or FM)
        if(type == 'ypr'){
            ##Plot
            ##op <- par(mfrow=c(1,1),new=F, mar = c(5, 4, 4, 4) + 0.3)

            ##plot Y/R & B/R
            label <- ifelse("Lc" %in% names(pes),  "Lc","tc")
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
            } ## only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area

            ## reference points
            if(length(N01) == 1){
                p.dat <- list_tc_Lc_runs[[tc_Lc_start]]
                py <- p.dat[,which(names(p.dat) == p.yield)]
                px <- p.dat[,which(names(p.dat) == p.FE)]
                ## F or  E 0.1
                if(N01 == Nmax){
                    add_shift <- 1.005
                }else  add_shift <- 1
                segments(x0 = -1, x1 = N01, y0 = py[which(px == N01)],
                         y1 = py[which(px == N01)],
                         col= 'darkgreen',lty = 1, lwd=1.5)
                segments(x0 = N01, x1 = N01, y0 = -1, y1 = py[which(px == N01)],
                         col= 'darkgreen',lty = 1, lwd=1.5)
                ## F or E max
                segments(x0 = -1, x1 = Nmax, y0 = py[which(px == Nmax)]*add_shift,
                         y1 = py[which(px == Nmax)]*add_shift,
                         col= 'goldenrod1',lty = 2, lwd=1.5)
                segments(x0 = Nmax*add_shift, x1 = Nmax*add_shift, y0 = -1, y1 = py[which(px == Nmax)],
                         col= 'goldenrod1',lty = 2, lwd=1.5)
                ## Legend
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


                ## current Exploitation rate or fishing mortality
                if(!is.null(pes$currents)){
                    currents <- pes$currents
                    if(!is.na(currents$curr.E)){
                        px <- ifelse(p.FE == "FM",currents$curr.F, currents$curr.E)
                        py <- ifelse(p.yield == "Y_R",currents$curr.YPR,currents$curr.YPR.rel)
                        points(px,py, pch = 16)
                        abline(v=px, col="grey30",lty=2)
                    }
                }
            }

            ##Biomass per recruit
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
                ## F or E 0.5
                segments(x0 = -1, x1 = N05, y0 = py2[which(px == N05)], y1 = py2[which(px == N05)],
                         col= 'red',lty = 3, lwd=1.5)
                segments(x0 = N05, x1 = N05, y0 = -1, y1 = py2[which(px == N05)],
                         col= 'red',lty = 3, lwd=1.5)
            }
            ##par(op)
        }

        ##Isopleths
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

                ## colours for plot
                pal <- colorRampPalette(rev(c(
                    rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1)
                )))
                if(yaxis1 == "B_R" | yaxis1 == "B_R.rel"){
                    pal <- colorRampPalette(c(
                        rgb(1,0.5,0.5), rgb(1,1,0.5), rgb(0.5,1,1), rgb(0.5,0.5,1)))
                }

                ##plot
                if(identify == TRUE){
                    dev.new(noRStudioGD = TRUE)
                }
                image(m, col=pal(100),
                      xlab = xlabel1, ylab = ylabel_iso)
                if(is.numeric(contour)){
                    contour(m, add=TRUE, nlevels = contour)
                }else if(contour == TRUE){
                    contour(m, add=TRUE)
                }

                ##mtext("Yield", line=0.5, side=3)
                if("currents" %in% names(pes) & mark){
                    points(x = curr_markX, y=curr_markY, pch=16, col="grey30")
                    abline(v=curr_markX, col="grey30",lty=2)
                    abline(h=curr_markY, col="grey30",lty=2)
                }
                if(identify == TRUE) image.identifier(m)
            }
        }
    }
}

#' @title Plotting production models
#'
#' @description This function plots CPUE and yield values against fishing effort
#'      resulting from the production models (\link{prod_mod}).
#'
#' @param x a object of the class \code{"prod_mod"}
#' @param display_MSY logical; should MSY be displayed in the graph?
#' @param ... optional parameters of plot function
#'
#' @examples
#' data(trawl_fishery_Java)
#' output <-  prod_mod(data = trawl_fishery_Java)
#' plot(output, display_years = TRUE)
#'
#' @importFrom graphics abline layout legend lines par plot points segments text title
#'
#' @references
#' Fox, W. W. Jr., 1970. An exponential surplus-yield model for optimizing exploited fish
#' populations. \emph{Trans.Am.Fish.Soc.}, 99:80-88
#'
#' Graham, M., 1935. Modern theory of exploiting a fishery and application to North Sea
#' trawling. \emph{J.Cons.CIEM}, 10(3):264-274
#'
#' Schaefer, M., 1954. Some aspects of the dynamics of populations important to the
#' management of the commercial marine fisheries. \emph{Bull.I-ATTC/Bol. CIAT}, 1(2):27-56
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @method plot prod_mod
#' @export

plot.prod_mod <- function(x, display_MSY = TRUE, ...){
    pes <- x
    f <- pes$f
    Y <- pes$Y
    fMSY.S <- pes$Schaefer_fMSY
    MSY.S <- pes$Schaefer_MSY
    fMSY.F <- pes$Fox_fMSY
    MSY.F <- pes$Fox_MSY
    a.F <- pes$Fox_lm[1]
    b.F <- pes$Fox_lm[2]
    a.S <- pes$Schaefer_lm[1]
    b.S <- pes$Schaefer_lm[2]
    CPUE.S <- pes$CPUE
    CPUE.F <- pes$ln_CPUE

    ##Plot
    op <- par(mfrow=c(2,1), xpd = FALSE,
              mar = c(1.2, 4, 1.5, 1) + 0.1,
              oma = c(5, 0.5, 2, 2) + 0.1)
    layout(matrix(c(1,2), nrow = 2, byrow=TRUE), heights=c(4,4))

    ## CPUE plot
    x <- seq(from = min(f),to = max(f),length.out = 500)
    y <- exp(a.F + b.F*x)
    plot(CPUE.S ~ f, xlab='', ylab='CPUE')
    lines(x = x, y = (a.S + x * b.S))
    ##abline(a = a.S, b=b.S)
    lines(x = x, y = y, col = 'blue', lty = 5)

    ## biomass plot
    x = seq(0,abs(a.S / b.S),length.out = 500)
    y.S = a.S*x + b.S*x^2
    y.F = x * exp(a.F + b.F*x)
    plot(x,rep(max(Y),length(x)),type='n',
         ylim=c(0,max(Y)*1.1), xlab = "", ylab = "Yield")
    points(Y ~ f)
    lines(x, y = y.S)
    lines(x, y = y.F,col='blue', lty = 5)
    if(display_MSY){
        ##segments(x0=fMSY.S,x1=fMSY.S,y0=0,y1=MSY.S,lty=2)
        ##segments(x0=-500,x1=fMSY.S,y0=MSY.S,y1=MSY.S,lty=2)
        points(fMSY.S,MSY.S,col='red', pch=18)
        text(fMSY.S, MSY.S,col='red',labels = expression(paste("MSY"[Schaefer])),
             adj = c(-0.05,1.1))
        ##segments(x0=fMSY.F,x1=fMSY.F,y0=0,y1=MSY.F,lty=2,col='blue')
        ##segments(x0=-500,x1=fMSY.F,y0=MSY.F,y1=MSY.F,lty=2,col='blue')
        points(fMSY.F,MSY.F, col='red', pch=18)
        text(fMSY.F, MSY.F, col='red',labels = expression(paste("MSY"[Fox])),
             adj = c(-0.05,1.1))
    }
    title(xlab = "Fishing effort", outer = TRUE, line = 2)
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend(x="top", col=c('black','blue'), bty ='n', horiz = TRUE,
           xpd = TRUE,legend = c('Schaefer','Fox'), lty=c(1,5))
    par(op)
}

#' @title Plotting time series production models
#'
#' @description This function plots objects of the class "prod_mod_ts".
#'
#' @param x a object of the class "prod_mod_ts",
#' @param correlation_plots logical; indicating if correlation plots should be displayed
#' @param ... additional parameters of the \link{plot} function
#'
#' @examples
#' data(emperor)
#' output <-  prod_mod_ts(emperor, method  = "Fox")
#' plot(output, correlation_plots = TRUE)
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics abline layout legend lines par plot
#' @importFrom stats lowess
#'
#' @references
#' Dharmendra, D., Solmundsson, J., 2005. Stock assessment of the offshore Mauritian banks
#' using dynamic biomass models and analysis of length frequency of the Sky Emperor
#' (\emph{Lethrinus mahsena}). Fisheries Training Program The United Nations University, 61
#'
#' @method plot prod_mod_ts
#' @export

plot.prod_mod_ts <- function(x, correlation_plots = FALSE, ...){
    pes <- x
    Y <- pes$Y
    yrs <- pes$year
    Bvec <- pes$Bvec
    CPUE <- pes$CPUE
    CPUE_hat <- pes$CPUE_hat
    K <- pes$K
    r <- pes$r
    q <- pes$q
    MSY <- pes$MSY
    Bmsy <- pes$Bmsy
    Emsy <- pes$Emsy
    method <- pes$method

    ## plotting
    if(!correlation_plots){
        op <- par(mfrow=c(2,2), xpd = FALSE,
                  mar = c(4, 4, 3, 1) + 0.1,
                  oma = c(4, 0.5, 1, 2) + 0.1)
        layout(matrix(c(1,2,3,4), nrow = 2, byrow=TRUE))

        ## Yield trajectory
        plot(yrs, Y, type="b", xlab="", ylab="Yield",
             main="Yield Trajectory",
             ylim=c(0, max(Y, na.rm = TRUE)*1.05))

        ## Biomass trajectory
        plot(yrs, Bvec, type="b", xlab="", ylab="Biomass",
             main="Biomass Trajectory", ylim=c(0, max(Bvec, na.rm = TRUE)*1.05))

        ## CPUE trajectory
        plot(yrs, CPUE, xlab="", ylab="CPUE",
             ylim=c(0, max(CPUE, na.rm = TRUE)*1.05), type="b",
             main = "CPUE Trajectory")
        lines(yrs, CPUE_hat, col=2, type="b")
        legend(x="bottomright",legend = c("observed", "predicted"),lty = 1, pch = 1,xpd = TRUE,
               col = c(1,2), cex = 0.8, bty = 'n', y.intersp = 0.8, x.intersp = 0.5)


### the equilibrium yield
        Blevels <- seq(0,ceiling(K),10)
        if(method == "Schaefer") EYlevels <- r*Blevels*(1-Blevels/K)
        if(method == "Fox") EYlevels <- r*Blevels*log(K/Blevels)
        EYlevels <- ifelse(EYlevels<0, 0, EYlevels)
        plot(Blevels, EYlevels, type="l" , xlab="Biomass",
             ylab="Yield", main="Yield Curve",
             xlim = c(0,max(c(Blevels,Bvec), na.rm = TRUE)),
             ylim = c(0,max(c(EYlevels,Y), na.rm = TRUE)))
        lines(Bvec,Y,type="b", col=2)
        abline(h=MSY,col=2)
        abline(v=Bmsy, col=2)

        par(op)
    }

    ## Validating Fit
    if(correlation_plots){
        ## plotting
        op <- par(mfrow=c(2,2), xpd = FALSE,
                  mar = c(4, 4, 4, 1) + 0.1,
                  oma = c(4, 0.5, 1, 2) + 0.1)
        layout(matrix(c(1,2,3,4), nrow = 2, byrow=TRUE))

        ## Correlation between yield and biomass
        plot(Y, Bvec, xlab="Yield", ylab="Biomass",
             main="Corr. between yield and biomass")
        lines(lowess(Y, Bvec), col=2)

        ## Correlation between biomass and CPUE
        plot(Bvec, CPUE, xlab="Biomass", ylab="CPUE",
             main="Corr. between biomass and CPUE")
        lines(lowess(Bvec, CPUE), col=2)

        ## Correlation between biomass and CPUE
        plot(log(Bvec), log(CPUE), xlab="Biomass", ylab="CPUE",
             main="Corr. between biomass and CPUE")
        lines(lowess(log(Bvec), log(CPUE)), col=2)

        ## observed vs predicted CPUE
        plot(CPUE, CPUE_hat, xlab="observed CPUE", ylab="predicted CPUE",
             main="Corr. observed and predicted CPUE")
        lines(lowess(CPUE, CPUE_hat), col=2)

        par(op)
    }
}


#' @title Millar's selectivity plot
                                        #
#' @description  This function plots the selectivity estimates of Millar's
#'    selectivity model (\code{\link{select_Millar}}).
#'
#' @param x a list of the class \code{"select_Millar"} containing the results of
#'    Millar's selectivity model
#' @param plotlens A vector with lengths which should be used for drawing the
#'    selection curves
#' @param standardise A parameter indicating if the retention should be realtive
#'    to the maximum value (Default: TRUE).
#' @param deviance_plot logical (Default: TRUE); indicating whether a plot of deviance residuals should
#'    be displayed
#' @param selectivity_plot logical (Default: TRUE); indicating whether a plot of relative retention selectivities
#' should be displayed
#' @param xlab_dev character string. Label for x axis of deviance plot. Default: "Length [cm]"
#' @param xlab_sel character string. Label for x axis of selectivity plot. Default: "Length [cm]"
#' @param ylab_dev character string. Label for y axis of deviance plot. Default: "Mesh size [cm]"
#' @param ylab_sel character string. Label for y axis of selectivity plot. Default: "Relative retention".
#' @param title_dev character string. Label for main title of deviance plot. Default: "Deviance residuals".
#' @param title_sel character string. Label for main title of selectivity plot. Default is
#' taken from the results of the select_Millar (e.g. res$rtype).
#' @param ... additional parameter options from plot function
#'
#' @examples
#' data(gillnet)
#'
#' output <- select_Millar(gillnet, x0 = c(60,4), rel.power = rep(1,8),
#'    rtype = "norm.loc", plot = FALSE)
#'
#' plot(output, plotlens = seq(40,90,0.1), deviance_plot = FALSE)
#'
#' @details This function draws a selectivity plot for the object class
#'    \code{"select_Millar"}, which is created by applying Millar's selectivity model
#'    \code{\link{select_Millar}}.
#'
#' @importFrom graphics abline axis matplot par plot points
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.
#'
#' @method plot select_Millar
#' @export


plot.select_Millar <- function(x,
                               plotlens = NULL,
                               standardise = TRUE,
                               deviance_plot = TRUE,
                               selectivity_plot = TRUE,
                               xlab_dev = "Length [cm]",
                               xlab_sel = "Length [cm]",
                               ylab_dev = "Mesh size [cm]",
                               ylab_sel = "Relative retention",
                               title_dev = "Deviance residuals",
                               title_sel = NULL,
                               ...
                               ){
    ## Adapted R code from Russell Millar (https://www.stat.auckland.ac.nz/~millar/selectware/)

    res <- x
    r <- rtypes_Millar(res$rtype)
    if(is.null(plotlens)) plotlens <- res$midLengths
    classes <- res$midLengths
    nlens <- length(classes)
    Dev.resids <- res$Dev.resids
    meshSizes <- res$meshSizes
    nmeshes <- length(meshSizes)
    AreLensUnique <- (length(classes)==length(unique(classes)))

    if(is.null(title_sel)){
        title_sel <- switch(res$rtype,
                            "norm.loc" = "Normal (common spread)",
                            "norm.sca" = "Normal",
                            "lognorm" = "Lognormal",
                            "binorm.sca" = "Bi-normal",
                            "bilognorm" = "Bi-lognormal",
                            "tt.logistic" = "Control and logistic",
                            ""
                            )
    }


    rmatrix <- outer(plotlens, res$meshSizes, r, res$par)
    rmatrix = t(t(rmatrix) * res$rel.power)

    if(standardise) rmatrix = rmatrix / max(rmatrix, na.rm = TRUE)

    ## define number of panels required for plots
    npanels <- sum(deviance_plot & ((nmeshes > 2 & AreLensUnique) | nmeshes == 2)) + sum(selectivity_plot)
    op <- par(mfrow = c(npanels,1))

    ## deviance plot
    if(deviance_plot & ((nmeshes > 2 & AreLensUnique) | nmeshes == 2)){

        if(nmeshes > 2 & AreLensUnique) {
            plot(1, 1, xlim=range(classes), xlab=xlab_dev, ylab=ylab_dev,
                 ylim=range(meshSizes)+(1/50)*c(-1,1)*(max(meshSizes)-min(meshSizes)), ## (cex/50)
                 yaxt="n", type="n", main=title_dev)
            axis(2,meshSizes,meshSizes,las=1)
            for(i in 1:nlens)
                for(j in 1:nmeshes)
                    points(classes[i],meshSizes[j],pch=ifelse(Dev.resids[i,j]>0,16,1),
                           cex=3*abs(Dev.resids[i,j])*1/(abs(max(Dev.resids))))   ## cex / (abs...)
        }else{
            if(nmeshes == 2) {
                Dev.resids.len=sign(Dev.resids[,2])*sqrt(apply(Dev.resids^2,1,sum))
                plot(classes, Dev.resids.len, type=ifelse(AreLensUnique, "h", "p"), las=1,
                     main=title_dev, xlab=xlab_dev, ylab=ylab_dev, cex=1)   ## cex = cex
                abline(h=0)
            }
        }

    }

    ## selectivity plot
    if(selectivity_plot){
        matplot(plotlens, rmatrix, type = "l", las = 1, ylim = c(0,1),
                xlab=xlab_sel, ylab=ylab_sel,
                main = title_sel, ...)
        abline(h = seq(0,1,0.25), lty = 3)
    }

    par(op)
}

#' @title Selectivity plot
                                        #
#' @description  This function plots the selectivity estimates of the
#'   function \code{\link{select}}.
#'
#' @param x a list of the class \code{"select"} containing the results of the
#'   gillnet selectivity function.
#' @param regression_fit logical; indicating if a plot with the fit of the regression
#'    line should be displayed
#' @param cols a specification for the two colours of the two selection curves.
#'    Default is c("darkgreen","orange").
#' @param ... additional parameters of the \link{plot} function
#'
#' @examples
#' data(tilapia)
#' output <- select(tilapia, plot = FALSE)
#' plot(output, regression_fit = TRUE)
#'
#' data(bream)
#' output <- select(bream, plot = FALSE)
#' plot(output, regression_fit = TRUE)
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics abline axis lines par plot segments text
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @method plot select
#' @export

plot.select <- function(x, regression_fit = TRUE,
                        cols = c("darkgreen","orange"), ...){

    if(length(cols) < 2) stop("Please provide two colours!")

    res <- x
    classes.num <- res$classes.num

    if(res$type == "gillnet"){
        numNet1 <- res$CatchPerNet_mat[,1]
        numNet2 <- res$CatchPerNet_mat[,2]
        msNet1 <- res$meshSizes[1]
        msNet2 <- res$meshSizes[2]
        LmNet1 <- res$LmNet1
        LmNet2 <- res$LmNet2
        s2 <- res$stand.dev
        SNet1 <- res$SNet1
        SNet2 <- res$SNet2
        linear_mod <- res$linear_mod
        lnNet2_Net1 <- res$lnNet2_Net1

        if(max(numNet2,na.rm=T) > max(numNet1,na.rm=T)){
            XX <- numNet2
        }else  XX <- numNet1

        classes.num.plot <- classes.num - 0.5
        xL <- seq(min(classes.num,na.rm=TRUE),max(classes.num,na.rm=TRUE),0.1)

        ##dev.new()
        ##create plot
        if(regression_fit){

            coeffs <- summary(linear_mod)$coefficients

            op <- par(mfrow = c(2,1), xpd = FALSE,
                      mar = c(4, 4, 2, 3) + 0.1,
                      oma = c(3, 0.5, 1, 2) + 0.1)

            plot(classes.num, lnNet2_Net1, xlab = "Length groups", ylab = "ln(Nnet2/Nnet1)",
                 ylim = c(min(lnNet2_Net1,na.rm=TRUE)*1.05,
                          max(lnNet2_Net1,na.rm=TRUE)*1.05))
            abline(a = coeffs[1,1], b = coeffs[2,1])

        }else op <- par(mfrow = c(1,1), mar = c(5, 5, 3, 3),
                        oma = c(3, 0.5, 1, 2))

        plot(classes.num.plot,XX,type='n', bty='n',xaxt ='n',
             ylab="Numbers caught",xlab="Fish length")
        lines(classes.num.plot,numNet1,type='s')
        axis(side=1,at=classes.num)
        lines(classes.num.plot,numNet2,type='s', lty=2)
        par(new=TRUE)
        plot(xL, exp(- ((xL - LmNet1)^2 / (2 * s2))), type = "l", col=cols[1],lwd=2.5,
             axes=F,bty = "n", xlab = "", ylab = "")
        lines(xL, exp(- ((xL - LmNet2)^2 / (2 * s2))), type = "l", col=cols[2],lwd=2.5,
              bty = "n", xlab = "", ylab = "")
        axis(side=4, at = pretty(range(SNet1)))
        mtext("Fractions retained", side=4, line=3)
        text(labels = paste("ms = ",msNet1,"cm"), x = LmNet1, y =
                                                                  ((exp(- ((LmNet1 - LmNet1)^2 / (2 * s2))))/0.92),
             col = cols[1],xpd = TRUE)
        text(labels = paste("ms = ",msNet2,"cm"), x = LmNet2, y =
                                                                  ((exp(- ((LmNet2 - LmNet2)^2 / (2 * s2))))/0.92),
             col = cols[2],xpd = TRUE)
        par(op)
    }

    if(res$type == "trawl_net"){
        reg.coeffs <- res$reg.coeffs
        S1 <- res$S1
        S2 <- res$S2
        SLobs <- res$SLobs
        L25 <- res$L25
        L50 <- res$L50
        L75 <- res$L75
        lnSL <- res$lnSL

        xL <- seq(min(classes.num,na.rm=TRUE),max(classes.num,na.rm=TRUE),0.1)

        ##dev.new()
        if(regression_fit){
            op <- par(mfrow = c(2,1), xpd = FALSE,
                      mar = c(4, 4, 2, 2) + 0.1,
                      oma = c(3, 0.5, 1, 1) + 0.1)
            plot(classes.num, lnSL, xlab = "Length groups", ylab = "ln(SL)",
                 ylim = c(min(lnSL,na.rm=TRUE)*1.05,
                          max(lnSL,na.rm=TRUE)*1.05))
            abline(a=reg.coeffs[1], b=reg.coeffs[2])

        }else op <- par(mfrow = c(1,1), mar = c(5, 5, 3, 3),
                        oma = c(3, 0.5, 1, 2))

        plot(SLobs ~ classes.num,
             xlab = "Fish length", ylab = "Fraction retained", pch = 16 , cex =1.4)
        lines(x = xL, y = 1/(1 + exp(S1 - S2 * xL)), col='blue', lwd = 1.5)
        segments(x0 = L25, x1 = L25, y0 = -1,
                 y1 = 1/(1 + exp(S1 - S2 * L25)),col= 'gray40',lty = 2, lwd=1.5)
        segments(x0 = 0, x1 = L25, y0 = 1/(1 + exp(S1 - S2 * L25)),
                 y1 = 1/(1 + exp(S1 - S2 * L25)),col= 'gray40',lty = 2, lwd=1.5)
        text(labels = "L25%", x = L25*1.042, y = (1/(1 + exp(S1 - S2 * L25)))/3,
             col = 'gray40')
        segments(x0 = L50, x1 = L50, y0 = -1,
                 y1 = 1/(1 + exp(S1 - S2 * L50)),col= 'gray40',lty = 2, lwd=1.5)
        segments(x0 = 0, x1 = L50, y0 = 1/(1 + exp(S1 - S2 * L50)),
                 y1 = 1/(1 + exp(S1 - S2 * L50)),col= 'gray40',lty = 2, lwd=1.5)
        text(labels = "L50%", x = L50*1.039, y = (1/(1 + exp(S1 - S2 * L50)))/3,
             col = 'gray40')
        segments(x0 = L75, x1 = L75, y0 = -1,
                 y1 = 1/(1 + exp(S1 - S2 * L75)),col= 'gray40',lty = 2, lwd=1.5)
        segments(x0 = 0, x1 = L75, y0 = 1/(1 + exp(S1 - S2 * L75)),
                 y1 = 1/(1 + exp(S1 - S2 * L75)),col= 'gray40',lty = 2, lwd=1.5)
        text(labels = "L75%", x = L75*1.035, y = (1/(1 + exp(S1 - S2 * L75)))/3,
             col = 'gray40')
        par(op)
    }
}


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
#' output <- VPA(hake, terminalE = 0.5, catch_unit = "'000")
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
#' @method plot VPA
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

    ##put together in dataframe
    if(yaxis == "numbers"){
        df.VPAnew <- data.frame(survivors = survivors,
                                nat.losses = natLoss,
                                catch = catch)
    }
    if(yaxis == "biomass" & dim(pes$plot_mat)[1] > 4){
        df.VPAnew <- data.frame(survivors = survivors * meanBodyWeight,
                                nat.losses = natLoss * meanBodyWeight,
                                catch = catch * meanBodyWeight)
    }else if(yaxis == "biomass" & dim(pes$plot_mat)[1] <= 4){
        stop("The information about the mean body weight per age or length class is missing!")
    }

    if(!display_last_class){
        df.VPAnew <- df.VPAnew[-dim(df.VPAnew)[1],]
        FM_calc <- FM_calc[-length(FM_calc)]
        classes.num <- classes.num[-length(classes.num)]
    }

    ##transpose matrix for barplot function
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

    ##save x axis positions
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

    ##create VPA plot
    ##dev.new()
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


