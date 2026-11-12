#' @title Univariate density plots of bootstrapping results
#'
#' @description Plots the univariate kernel density estimates of the VBGF
#'   parameters from bootstrapping results (\link{ELEFAN_GA_boot},
#'   \link{ELEFAN_SA_boot}) together with the confidence interval and the
#'   parameter value with maximum density.
#'
#' @param res list of class `lfqBoot`
#' @param CI numeric. Confidence interval level (Default: 95)
#' @param use_hist logical. Plot histogram in addition to smoothed kernel density.
#' @param nbreaks numeric. Number of breaks in the histogram.
#' @param truePar if true parameters known lines can be added to graph
#' @param display_val logical. Should the values of the CI limits and maximum
#'   density be displayed (Default: `TRUE`)?
#' @param display_median logical. Should the median be displayed (Default: `FALSE`)?
#' @param display_legend logical. Should a legend be displayed (Default: `FALSE`)?
#' @param display_rugs logical. Should the individual estimates be displayed as
#'   rugs (Default: `FALSE`)?
#' @param mar vector. Inner margins settings. See `?par`.
#' @param oma vector. Outer margins settings.See `?par`.
#' @param mgp vector. See `?par`.
#' @param tcl vector. See `?par`.
#' @param cex numeric. See `?par`.
#' @param mfrow vector. See `?par` (default: NA).
#' @param ylim optional list with y-axis limits, one vector per parameter
#'   (Default: NA).
#' @param ... Additional arguments passed to `par`.
#'
#' @return plot
#'
#' @importFrom graphics par plot hist rect polygon lines segments abline text
#' @importFrom graphics mtext box plot.new legend
#' @importFrom grDevices rgb
#' @importFrom stats na.omit median
#'
#' @examples
#' data(alba_boot)
#' univariate_density(alba_boot)
#'
#' @export
#'
univariate_density <- function(res, CI=95, use_hist = FALSE, nbreaks = 10,
                               truePar = NA, display_val=TRUE, display_median=FALSE,
                               display_legend = FALSE, display_rugs = FALSE,
                               mar = c(1.5,2,2,0), oma = c(1.5,0,0,0.5),
                               mgp = c(2,0.5,0), tcl = -0.25, cex = 1, mfrow=NA,
                               ylim=NA,...
                               ){
    res <- res$bootRaw
    na <- ifelse(display_legend,ncol(res)+1,ncol(res))
    nc <- min(c(6,na))
    nr <- ceiling(na/nc)
    if(is.na(mfrow[1])) mfrow <- c(nr,nc)
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    par(
        mfcol = c(1, ncol(res)),
        mar = mar, oma = oma, mfrow = mfrow,
        mgp = mgp, tcl = tcl, cex = cex, ...
    )

    VARS <- list(
        Linf = expression(bolditalic(L)[infinity]),
        K = expression(bolditalic(K)),
        t_anchor = expression(bolditalic(t)[anchor]),
        C = expression(bolditalic(C)),
        ts = expression(bolditalic(t)[s]),
        phiL = expression(bolditalic("phi'")),
        agemax = expression(bolditalic(t)[max]),
        M_Then = expression(bolditalic(M)),
        M_Pauly = expression(bolditalic(M)),
        Z = expression(bolditalic(Z)),
        E = expression(bolditalic(E)),
        FM = expression(bolditalic(F)),
        L50 = expression(bolditalic("L50")),
        L75 = expression(bolditalic("L75")),
        N = expression(bolditalic(N)),
        B = expression(bolditalic(B)),
        F01 = expression(bolditalic("F01")),
        Fmax = expression(bolditalic(Fmax)),
        F05 = expression(bolditalic("F05")),
        FF01 = expression(bolditalic("FF01")),
        FFmax = expression(bolditalic(FFmax)),
        FF05 = expression(bolditalic("FF05")),
        ypr = expression(bolditalic(ypr)),
        yprrel = expression(bolditalic(yprrel)),
        bpr = expression(bolditalic(bpr)),
        bprrel = expression(bolditalic(bprrel))
    )


    ## univariate plots
    for(i in seq(ncol(res))){
        tmp <- as.numeric(na.omit(res[,i]))
        x <- try(ks::kde(tmp),silent=TRUE)
        if(inherits(x, "try-error")){
            ## e.g. constant values (no density can be estimated)
            plot(seq(0,1,0.2),
                 seq(min(tmp)*0.9,max(tmp)*1.1,length.out = 6),
                 t="n",
                 xlim = c(0,1),
                 xaxs = "i",
                 ylab="", xlab="", col=1, lty=1
                 )
        }else{
            h = hist(tmp, plot=FALSE, breaks = nbreaks)
            xlim <- c(0, max(x$estimate))
            if(use_hist){
                xlim <- range(c(xlim, max(h$density)))
            }
            xlim <- xlim * c(0,1.1)

            ylims <- range(x$eval.points)
            if(any(!is.null(ylim)) && any(!is.na(ylim))){
                if(!is.na(ylim[[i]][1])){
                    ylims[1] <- ylim[[i]][1]
                    ylims[2] <- ylim[[i]][2]
                }
            }


            plot(x$estimate, x$eval.points, t="n",
                 xlim = xlim,
                 ylim = ylims,
                 xaxs = "i",
                 ylab="", xlab="", col=1, lty=1
                 )
            usr <- par()$usr

            CItxt <- paste0(round(100-CI), "%")
            inCI <- rle( x$estimate > x$cont[CItxt] )
            start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
            end.idx <- cumsum(inCI$lengths)
            limCI <- range(x$eval.points[start.idx[min(which(inCI$values))]:end.idx[max(which(inCI$values))]])

            in1 <- which(x$estimate > x$cont["99%"])
            mean1 <- mean(x$eval.points[in1])

            if(use_hist){
                rect(
                    xleft = 0, ybottom = h$breaks[-length(h$breaks)],
                    xright = h$density, ytop = h$breaks[-1],
                    col = "grey90", border = "grey50"
                )
            }else{
                for(j in seq(inCI$lengths)){
                    if(inCI$values[j]){
                        polygon(
                            y = c(x$eval.points[start.idx[j]:end.idx[j]], rev(x$eval.points[start.idx[j]:end.idx[j]])),
                            x = c(x$estimate[start.idx[j]:end.idx[j]], rep(0, length(x$estimate[start.idx[j]:end.idx[j]]))),
                            col = "grey90",
                            border = NA, lty = 3, lwd = 1
                        )
                    }
                }
            }
        }

        if(!inherits(x, "try-error")){
            lines(x$estimate, x$eval.points, lwd = 1, col = "grey50")

            ## rug
            if(display_rugs)
                segments(x0 = 0, x1 = par()$cxy[1]*0.3, y0 = x$x, y1 = x$x, col=rgb(0,0,0,0.5), lwd=0.3)

            ## true pars
            if(!is.na(truePar[1])){
                trupari <- truePar[match(names(res)[i], names(truePar))]
                abline(h=trupari, lty=2, col=2)
            }

            ## range of CI
            abline(h = limCI, lty = 2, lwd=1, col = 1)
            if(display_val) text(y =c(limCI), x = mean(usr[1:2]),
                                 labels = paste(sprintf("%.2f", round(c(limCI),2))),
                                 pos = c(1,3), offset = 0.25, col = 1)

            ## median
            if(display_median) abline(h=median(x$x), col=4, lty=2)

            ## maximum density
            abline(h = mean1, lty = 1, lwd=1, col = 1)
            if(display_val) text(y =  mean1, x = mean(usr[1:2]),
                                 labels = sprintf("%.2f", round(mean1,2)),
                                 pos = 3,
                                 offset = 0.25, col = 1
                                 )


        }else{
            abline(h=unique(tmp), lty = 3)
        }
        varlab <- VARS[[match(names(res)[i], names(VARS))]]
        if(is.null(varlab)) varlab <- names(res)[i]
        mtext(varlab, line=0.25, side=3, font = 2)
        box(lwd=1.2)
    }
    mtext("Density", side = 1, line = 0, outer = TRUE)

    if(display_legend){
        plot.new()
        legend("center",
               legend=c("true","median","max density","CI"),
               lty=c(2,2,1,2), col=c(2,4,1,1), bty="n")
    }
}
