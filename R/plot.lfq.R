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
#'   \item \strong{Linf} length infinity,
#'   \item \strong{K} curving coefficient,
#'   \item \strong{C} amplitude of growth oscillation (optional),
#'   \item \strong{WP} winter point (WP = ts + 0.5) (optional);
#' }
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; should positive peaks be flagged out?
#'    (Default : TRUE)
#' @param col colour of growth curves (default: 1)
#' @param add.image logical; should and image be added to the lfq plot
#'    (default : TRUE)
#' @param col.image colour of image, by default (NULL) red and blue colours
#'    are used
#' @param zlim the minimum and maximum z values for which colors should be
#'    plotted (default : NULL).
#' @param zlimtype indicating if zlim should be based on the range of
#'    the catch ("range") or based on the maximum catch ("balanced", default). This parameter is only
#'    considered if zlim is NULL.
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
#' data(synLFQ4)
#' res <- lfqRestructure(synLFQ4)
#' plot(x = res, Fname = "rcounts")
#' plot(res, Fname = "rcounts", par = list(Linf = 80, K = 0.5, t_anchor = 0.25))
#' plot(res, Fname = "catch", par = list(Linf = 80, K = 0.5,
#'    t_anchor = 0.25, C= 0.75, WP = -0.5))
#'
#' # plot length frequency data without restructuring
#' class(synLFQ4) <- "lfq"
#' plot(synLFQ4, Fname = "catch")
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
#' @export

plot.lfq <- function(x, Fname = "rcounts",  # alternative : "catch"
                     par = NULL,
                     agemax = NULL,
                     flagging.out = TRUE,
                     col = "blue",
                     add.image = TRUE,
                     col.image = NULL,
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

  bin.width <- diff(classes)
  bin.lower <- classes - c(bin.width[1], bin.width)/2
  bin.upper <- classes + c(bin.width, bin.width[length(bin.width)])/2

  # bin height scaling
  if(Fname == "catch"){
    sc <- min(diff(dates)) * 0.8 / max(abs(catch))
  }
  if(Fname == "rcounts"){
    sc <- min(diff(dates)) * 0.5 / max(abs(catch))
  }


  # image colour
  if(is.null(col.image)){
    pal <- colorRampPalette(c(rgb(1,0.8,0.8), rgb(1,1,1), rgb(0.8,0.8,1)))
    col.image <- pal(21)
  }

  # zlim value + type
  if(is.null(zlim) & zlimtype == "balanced"){
    zlim = c(-1,1) * max(abs(catch), na.rm=TRUE)
  }
  if(is.null(zlim) & zlimtype == "range"){
    zlim = range(catch, na.rm = TRUE)
  }

  # with image
  if(add.image){
    #par(mar = c(5, 5, 1, 1) + .1)
    image(x=dates, y=classes, z=t(catch), col=col.image, zlim=zlim,
          xaxt="n",  xlab = xlab, ylab = ylab, ...)
  }
  # without image
  if(add.image == FALSE){
    image(x=dates, y=classes, z=t(catch), col=NA, zlim=zlim,
          xaxt="n",  xlab = xlab, ylab = ylab, ...)
  }

  # add time axis
  if(date.axis == "modern"){
    axis.Date(side = 1, x=dates, at=date.at, format = date.format)
  }else if(date.axis == "traditional"){
    axis.Date(side = 1, x = dates, at = date.at, format = "%b")
    year <- unique(format(dates, "%Y"))
    date_seq <- seq.Date(dates[1],dates[length(dates)], by = "month")
    date_label <- format(date_seq, "%m")
    year_pre <- which(date_label %in% "01")
    if(!(1 %in% year_pre)) year_pre <- c(1,which(date_label %in% "01"))
    dates_for_years <- as.Date(paste(format(date_seq,"%Y"),date_label,"01",sep="-"))
    year_ticks <- dates_for_years[year_pre]
    mtext(side = 1, at = year_ticks, text = year, line = 2.5)
  }

  # frame
  box(col = "gray40") #bty = "L"

  #Histograms
  for(i in seq(length(dates))){
    score.sc <- catch[,i] * sc
    for(j in seq(classes)){
      polygon(x = c(dates[i], dates[i], dates[i]-score.sc[j], dates[i]-score.sc[j]),
              y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
              col = c("white", "black")[(score.sc[j]>0)+1],
              border = "grey20", lwd = 1)
    }
  }

  if("par" %in% names(x) & is.null(par) & draw){
    Lt <- lfqFitCurves(lfq = x, par = x$par,
                       agemax = x$agemax, flagging.out = flagging.out, draw = TRUE)
  }
  if(!is.null(par) & draw){
    Lt <- lfqFitCurves(x, par = par,
                       agemax = agemax, flagging.out = flagging.out, draw = TRUE)
  }
}
