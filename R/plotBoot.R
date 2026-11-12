#' @title Plot length-frequency data with bootstrapped growth curves
#'
#' @description Plots the length-frequency data together with the growth curves
#'   with maximum density and their confidence interval from bootstrapping
#'   results (\link{ELEFAN_GA_boot}, \link{ELEFAN_SA_boot}).
#'
#' @param lfq a length frequency object of the class `lfq`
#' @param boot an object of class `lfqBoot`
#' @param Fname name of the element of `lfq` to plot, either `"catch"`
#'   (default) or `"rcounts"` (restructured data, see \link{lfqRestructure})
#' @param agemax numeric. Maximum age of the growth curves. By default, derived
#'   from the bootstrapped estimates of K.
#' @param rel logical; should relative frequencies (per sample) be displayed?
#'   (Default: `FALSE`)
#' @param hist.sc numeric; scaling of the histograms (Default: 0.5)
#' @param hist.col colours of the histograms for negative and positive values
#' @param image.col colours of the background image
#' @param region.col optional; colour of the plot region (replaces the
#'   background image)
#' @param zlim limits of the background image
#' @param zlimtype either `"balanced"` (default) or `"range"`
#' @param date.axis either `"traditional"` (default) or `"modern"`
#' @param date.at dates for the date axis
#' @param date.format format of the dates for `date.axis = "modern"`
#' @param xlab label for x-axis
#' @param ylab label for y-axis (Default: `"Length classes"`)
#' @param CI numeric; confidence interval level (Default: 95)
#' @param ci.col colour of the confidence interval lines
#' @param ci.lty line type of the confidence interval lines
#' @param ci.lwd line width of the confidence interval lines
#' @param maxd.col colour of the growth curves with maximum density
#' @param maxd.lty line type of the growth curves with maximum density
#' @param maxd.lwd line width of the growth curves with maximum density
#' @param ... additional arguments passed to \code{\link[graphics]{image}}
#'
#' @return plot
#'
#' @importFrom graphics image rect axis.Date mtext polygon lines box par
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom stats cov quantile
#'
#' @examples
#' data(alba)
#' data(alba_boot)
#' plotBoot(lfq = alba, boot = alba_boot)
#'
#' # restructured data
#' plotBoot(lfq = lfqRestructure(alba, MA = 7), boot = alba_boot, Fname = "rcounts")
#'
#' @export
#'
plotBoot <- function(lfq,
                     boot,
                     Fname = "catch",
                     agemax = NULL,
                     rel = FALSE,
                     hist.sc = 0.5,
                     hist.col = c("white", "black", "orange","darkgreen"),
                     image.col = NULL,
                     region.col = NULL,
                     zlim = NULL,
                     zlimtype = "balanced",   # alternative : "range"
                     date.axis = "traditional",  # alternative : "modern"
                     date.at = seq(as.Date("1500-01-01"), as.Date("2500-01-01"), by="months"),
                     date.format = "'%y-%b", xlab = "", ylab = "Length classes",
                     CI = 95,
                     ci.col = 1, ci.lty = 2, ci.lwd = 1,
                     maxd.col = 1, maxd.lty = 1, maxd.lwd = 2,
                     ...){

    dates <- lfq$dates
    classes <- lfq$midLengths
    catch <- get(Fname, lfq)

    ## display relative catches (relative to number of samples per month)
    if(rel){
        catchrel <- catch
        for(i in 1:ncol(catch)){
            catchrel[,i] <- catch[,i]/colSums(catch, na.rm = TRUE)[i]
        }
        catch <- catchrel
        catch[is.nan(catch)] <- 0
    }
    ## bin height scaling
    sc <- unclass(min(diff(dates)) * hist.sc / max(abs(catch)))

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
    image(
        x=dates, y=classes, z=t(catch), col=image.col, zlim=zlim,
        xaxt="n", xlab = xlab, ylab = ylab, ...)

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
    for(i in seq(length(dates))){
        score.sc <- catch[,i] * sc
        for(j in seq(classes)){
            polygon(
                x = c(dates[i], dates[i], dates[i]-score.sc[j], dates[i]-score.sc[j]),
                y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
                col = hist.col[(score.sc[j]>0)+1],
                border = "grey20", lwd = 1)
        }
    }

    ## plotting cohort growth curves with uncertainty
    res <- boot$bootRaw
    if(is.null(agemax))
        agemax <- max(ceiling((1/-res$K)*log(1-((res$Linf*0.95)/res$Linf))))

    ## VBGF parameters only
    colind <- which(names(res) %in% c("Linf","K","t_anchor"))
    if("C" %in% names(res)) colind <- c(colind, which(names(res) %in% c("C","ts")))
    x <- res[,colind]

    ## First a fitting round to look for shifts in t_anchor
    age = seq(0, agemax, 0.01)
    Lt0 <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    Lt_minus <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    Lt_plus <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    for(i in seq(ncol(Lt0))){
        par0 <- par_minus <- par_plus <- as.list(x[i,])
        par_minus$t_anchor <- par0$t_anchor-1
        par_plus$t_anchor <- par0$t_anchor+1
        Lt0[,i] <- VBGF(param = par0, t = age)
        Lt_minus[,i] <- VBGF(param = par_minus, t = age)
        Lt_plus[,i] <- VBGF(param = par_plus, t = age)
    }

    ## replace negative lengths with NA
    Lt0 <- replace(Lt0, Lt0<0, NaN)
    Lt_minus <- replace(Lt_minus, Lt_minus<0, NaN)
    Lt_plus <- replace(Lt_plus, Lt_plus<0, NaN)

    ## determine if a positive or negative shift improves overall covariance for each permutation
    cov0 <- cov(Lt0, use = "pair")
    shift <- 0 * seq(ncol(Lt0))
    for(i in seq(ncol(Lt0))){
        cov_minus <- cov(x = Lt_minus[,i], y = Lt0, use = "pair")
        cov_plus <- cov(x = Lt_plus[,i], y = Lt0, use = "pair")
        shift[i] <- c(0,-1,1)[which.max(c(sum(cov0[i,]), sum(cov_minus), sum(cov_plus)))]
    }

    ## fixed predictions
    agenew <- seq(0+min(shift), agemax+max(shift), 0.01)
    Lt <- matrix(NaN, nrow=length(agenew), ncol=nrow(x))
    for(i in seq(ncol(Lt))){
        par0 <- as.list(x[i,])
        par0$t_anchor <- par0$t_anchor + shift[i]
        Lt[,i] <- VBGF(param = par0, t = agenew)
    }

    ## multivariate kernel density estimate for max. density
    H <- ks::Hpi(x, nstage = 1)
    fhat <- ks::kde(x = x, H = H, eval.points = x)

    ## maximum density
    max_dens <- as.list(fhat$eval.points[which.max(fhat$estimate),])

    ## determine which resamples are in the CI
    inCI <- fhat$estimate > quantile(fhat$estimate, probs = (100-CI[1])/100 )
    limCI <- data.frame(t = agenew, min = NaN, max = NaN)
    for(i in seq(agenew)){
        limCI$min[i] <- min(Lt[i, which(inCI) ], na.rm = TRUE)
        limCI$max[i] <- max(Lt[i, which(inCI) ], na.rm = TRUE)
    }

    ## one cohort per year; VBGF with t_anchor already starts the curves at
    ## t_anchor, so that the cohort of a year starts at year + t_anchor (as in
    ## lfqFitCurves)
    years <- as.numeric(as.character(unique(format(dates, "%Y"))))
    years2 <- seq(min(years) - agemax, max(years), 1)
    for(i in seq_along(years2)){
        xvals <- yeardec2date(years2[i] + agenew)

        lines(xvals, limCI$max, col = ci.col, lwd = ci.lwd, lty = ci.lty)
        lines(xvals, limCI$min, col = ci.col, lwd = ci.lwd, lty = ci.lty)

        lines(xvals, VBGF(param = max_dens, t = agenew),
              col = maxd.col, lwd = maxd.lwd, lty = maxd.lty)
    }

    ## frame
    box(lwd=1.2)
}
