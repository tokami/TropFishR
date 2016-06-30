#' @title Plotting of length frequency data (with VBGF curves)
#'
#' @description This function plots length frequency (lfq) samples sequentially arranged in time. An object
#'  of "lfq" class is obatined by applying the \code{\link{lfqRestructure}} function. In case growth
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
#' @param multiple_best_fits numeric; indicating which startingPoints should be used in case there are more
#'    than one (default: 1)
#' @param col colour of growth curves (default: "blue")
#' @param ... additional options of the plot function
#'
#'
#' @examples
#' data(trout)
#' res <- lfqRestructure(trout)
#' plot(x = res, Fname = "rcounts")
#' plot(res, Fname = "rcounts", par = list(Linf = 16, K = 0.77))
#' plot(res, Fname = "catch", par = list(Linf = 16, K = 0.77, C= 0.5, WP = 0.6))
#'
#' # plot length frequency data without restructuring
#' class(trout) <- "lfq"
#' plot(trout, Fname = "catch")
#'
#' @details This function runs \code{\link{lfqFitCurves}} when growth parameters are provided. Thereby the starting point
#'    calculated. It can happen that different starting points return the same fit (ESP value), with
#'    \code{multiple_best_fits} the preferred starting point can be chosen.
#'
#' @import lubridate
#' @importFrom grDevices colorRampPalette dev.new dev.off recordPlot rgb
#' @importFrom graphics abline axis box lines mtext par plot rect
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

plot.lfq <- function(x, Fname = "rcounts",
                     par = NULL,
                     multiple_best_fits = 1, col = "blue", ...){
  dates <- x$dates
  classes <- x$midLengths
  catch <- get(Fname, x)
  n_samples <- dim(catch)[2]

  interval <- (classes[2]-classes[1])/2

  # get length of smapling period
  # continuous time in years
  julian_days1 <- as.numeric(format(dates, format="%Y")) + as.numeric(format(dates, format="%j"))/366
  days.years <- julian_days1 - julian_days1[1]  # OLD: #days <- as.numeric(dates - as.Date((dates[1]))) #days.years <- days/365

  # for x axis
  months_first <- update(object = dates, mdays = 1, hours = 0, minutes = 0, seconds = 0)
  date_seq <- seq.Date(months_first[1],months_first[length(months_first)], by = "month")
  date_label <- format(date_seq, "%b")
  julian_days2 <- as.numeric(format(date_seq, format="%Y")) + as.numeric(format(date_seq, format="%j"))/366
  date_ticks <- julian_days2 - julian_days1[1]  #date_ticks <- as.numeric(date_seq - as.Date((dates[1])))/365

  # for years beneath x axis
  year <- unique(format(dates, "%Y"))
  year_num <- floor(date_ticks[length(date_ticks)] - date_ticks[1]) + 1
  year_pre <- which(date_label %in% "Jan")
  if(!(1 %in% year_pre)) year_pre <- c(1,which(date_label %in% "Jan"))
  #year_pre <- seq(1,1:year_num*12,length.out = year_num)
  year_ticks <- date_ticks[year_pre]

  #for plotting
  xplot <- c((min(days.years, na.rm = TRUE)-0.2),(max(days.years, na.rm = TRUE)+0.2))
  yplot <- c(0,classes[length(classes)]*1.1)  #### c(classes[1],classes[length(classes)])


  # Plot with rearranged histogramms
  par(mar = c(5, 5, 1, 1) + .1)
  plot(xplot, yplot, type = "n",axes = FALSE, #ylim = c(0,max(yplot)),
       ann = FALSE,  xaxs = "i", yaxs = "i")
  axis(1,at = date_ticks, labels = date_label, cex.axis = 1.2)
  mtext(side = 1, at = year_ticks, text = year, line = 3, cex = 1.2)
  axis(2, cex.axis = 1.2)
  mtext(side = 2, outer = F, line = 3.5, "Length [cm]", cex = 1.2)
  box(col = "gray40") #bty = "L"

  min_delta_t <- min(diff(days.years),na.rm=TRUE)
  max_peak <- max(abs(catch),na.rm = TRUE)
  if(min_delta_t <= max_peak){
    bar_height_adjust <- floor(max_peak /  min_delta_t)
  }else bar_height_adjust = 1


  #Histograms
  y.b <- classes - interval
  y.t <- classes + interval
  for(histi in 1:n_samples){
    x.l <- rep(days.years[histi],length(classes))
    x.r <- x.l - catch[,histi] / bar_height_adjust
    abline(v = x.l[1],col = "gray70")
    # positive ones
    posi <- x.r > x.l
    rect(xleft = x.l[posi], ybottom = y.b[posi], xright = x.r[posi], ytop = y.t[posi],
         col = "white", border = "black")
    # negative ones
    negi <- x.r <= x.l
    rect(xleft = x.l[negi], ybottom = y.b[negi], xright = x.r[negi], ytop = y.t[negi],
         col = "black", border = "black")
  }

  # when parameters given then plot growth curves
  if(!is.null(par[[1]]) | ("K_opt" %in% names(x) & "Linf_fix" %in% names(x))){
    if(!is.null(par[[1]])){
      if(!("Linf" %in% names(par)) | !("K" %in% names(par))) stop("At least 'Linf' and 'K' have to be defined in par!")
      Linfi <- get("Linf",par)
      Ki <- get("K",par)
      if("C" %in% names(par)){
        C <- get("C",par)
      }else C <- 0
      if("WP" %in% names(par)){
        WP <- get("WP",par)
      }else WP <- 0
      if("ts" %in% names(par)) WP <- 0.5 + get("ts",par)
    }
    if(("K_opt" %in% names(x) & "Linf_fix" %in% names(x))){
      Linfi <- x$Linf_fix
      Ki <- x$K_opt[multiple_best_fits]
      C <- x$C
      WP <-  x$WP
      par <- list(Linf = Linfi, K = Ki, C = C, WP = WP)
    }

    res <- lfqFitCurves(lfq = x, par = par)

    if(length(res$startingSample) == 1){
      startingSample <- res$startingSample
      rel_time_startingSample <- days.years[as.numeric(startingSample)]
      startingLength <- as.numeric(res$startingLength)
    }else{
      startingSample <- res$startingSample[multiple_best_fits]
      rel_time_startingSample <- days.years[as.numeric(startingSample)]
      startingLength <- as.numeric(res$startingLength[multiple_best_fits])
    }

    # maximum age to be displayed - end of growth curve
    Lmax <- Linfi * 0.98
    tmax_plot <- log(1-Lmax/Linfi) / -Ki

    # lookup table for soVBGF
    if(C != 0 | WP != 0.5){
      lookup_age <- seq(0,tmax_plot+30,0.01)   ##ISSUE: tmax
      lookup_length <- Linfi * (1 - exp(-Ki * lookup_age + (((C*Ki)/(2*pi)) * sin(2*pi*(lookup_age-(WP-0.5)))) - (((C*Ki)/(2*pi))*sin(2*pi*(-(WP-0.5))))))
    }

    tsample_t0 <- log(1-startingLength/Linfi) / -Ki
    if(C != 0 | WP != 0.5){
      lookup_ind <- which.min(abs(lookup_length - startingLength))
      tsample_t0 <- lookup_age[lookup_ind]
    }
    tstart <- -tsample_t0 + rel_time_startingSample

    # maximum age to be displayed - end of growth curve - get more exact value if seasonalised
    if(C != 0 | WP != 0.5){
      lookup_ind <- which.min(abs(lookup_length - Lmax))
      tmax_plot <- lookup_age[lookup_ind]
    }


    #number of cohorts
    n_cohorti <- ceiling(abs(tmax_plot))
    n_cohorti <- ifelse(n_cohorti == 0,1,n_cohorti)
    n_cohorti <- ifelse(is.na(n_cohorti),1,n_cohorti)


    # x axis intercept of growth curves in comparison to:  days.years to get this sequence
    # negative tstart is intercept with x axis of starting cohort
    younger_cohorts <- floor(abs(tstart - days.years[length(days.years)]))
    older_cohorts <- floor(tmax_plot - tstart)

    add_t_plot_seq <- seq(-older_cohorts,younger_cohorts,1)

    births_cohorts <- tstart + add_t_plot_seq
    deaths_cohorts <- births_cohorts + tmax_plot

    for(cohorti in 1:length(add_t_plot_seq)){
      addi <- add_t_plot_seq[cohorti]
      t_cohorti <- seq(births_cohorts[cohorti],deaths_cohorts[cohorti],0.005)
      Lt <- Linfi * (1 - exp(-Ki * (t_cohorti - tstart - addi) +
                               (((C*Ki)/(2*pi)) * sin(2*pi*((t_cohorti - tstart - addi)-(WP-0.5)))) -
                               (((C*Ki)/(2*pi)) * sin(2*pi*(0-(WP-0.5))))))
      lines(y = Lt, x = t_cohorti, lty=2, col=col)
    }
    return(res)
  }
}

