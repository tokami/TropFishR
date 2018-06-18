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
#'    displayed (default: "catch")
#' @param par a list with following growth parameters (default NULL):
#'  \itemize{
#'   \item \strong{Linf} asymptotic length,
#'   \item \strong{K} growth coefficient,
#'   \item \strong{t_anchor} time at length zero,
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
#' plot(res, Fname = "rcounts", par = list(Linf = 14, K = 1.1, t_anchor = 0.3))
#'
#' # add soVBGF curves, adjust hist.sc and xlim
#' plot(res, Fname = "catch", curve.col=4,
#'   par = list(Linf = 14, K = 1.1, t_anchor = 0.3, C = 0.2, ts = 0.75),
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
  Fname = "catch",  # alternative : "rcounts"
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
  zlimtype = "balanced",   # alternative : "range"
  date.axis = "traditional",  # alternative : "modern"
  date.at = seq(as.Date("1500-01-01"), as.Date("2500-01-01"), by="months"),
  date.format = "'%y-%b", xlab = "", ylab = "Length classes",
  draw = TRUE,
  ...){

    dates <- x$dates
    classes <- x$midLengths
    catch <- get(Fname, x)
    agemax <- ifelse("agemax" %in% names(x), x$agemax, NA)
    spawningTimes <- ifelse("spawningTimes" %in% names(x), x$spawningTimes, NA)        
    if("par" %in% names(x)){
        if(length(x$par$t_anchor) == 2){
            spawningTimes <- 2
        }
    }


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

    # image colour
    if(is.null(image.col)){
      pal <- colorRampPalette(c(rgb(1,0.8,0.8), rgb(1,1,1), rgb(0.8,0.8,1)))
      image.col <- pal(21)
    }
    if(!is.null(region.col)){
      image.col <- NA
    }

    # zlim value + type
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

    # add time axis
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
        # if(score.sc[j] != 0){
          polygon(
            x = c(dates[i], dates[i], dates[i]-score.sc[j], dates[i]-score.sc[j]),
            y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
            col = hist.col[(score.sc[j]>0)+1],
            border = "grey20", lwd = 1)
        # }
      }
    }
    }

    # optional addition of cohort growth curves
    if("par" %in% names(x) & is.null(par) & draw){
      Lt <- lfqFitCurves(lfq = x, par = x$par,
                         agemax = agemax, draw = TRUE, col=curve.col,
                         spawningTimes = spawningTimes)
    }
    if(!is.null(par) & draw){
        Lt <- lfqFitCurves(x, par = par,
                           agemax = agemax, draw = TRUE, col=curve.col,
                           spawningTimes = spawningTimes)
    }

    # frame
    box()
}
