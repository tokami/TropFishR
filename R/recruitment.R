#' calculate time for new length
#'
#' @param par VBGF parameters. If soVBGF, then ts should relate
#' to time of calendar year rather than time of age year.
#' @param Ltnew New length for which time is calculated
#' @param Lt Known length at time t
#' @param t known time t
#' @param tincr time increment (year fraction) used for iterative
#' hind-/fore-casting. Use smaller values for increased precision
#' (Default: tincr = 0.01)
#'
#' @return numeric value representing year decimal
#' @export
#'
#' @examples
#'
#' data(synLFQ4)
#' synLFQ4$par <- list(Linf = 80, K = 0.5, ta = 0.25, C = 0.75, ts = 0.5)
#'
#' # time at length = 0
#' calc_tnew(par = synLFQ4$par, Lt = 70, t = 2000, Ltnew = 0)
#'
#' # time at length = 75
#' calc_tnew(par = synLFQ4$par, Lt = 70, t = 2000, Ltnew = 75)
#'
calc_tnew <- function(
  par = NULL,
  Ltnew = NULL,
  Lt = NULL,
  t = NULL,
  tincr = 0.01
){
  if(Ltnew == Lt | Ltnew >= par$Linf | Lt >= par$Linf | Lt < 0){
    tnew <- NA
  }else{
    Lt.i <- Lt
    if(class(t)=="Date"){
      t.i <- date2yeardec(t)
    }else{
      if(class(t)=="numeric"){
        t.i <- t
      }else{
        stop("Error: t must be of class 'Date' or 'numeric'" )
      }
    }
    seasonalized <- !is.null(par$C)
    if(!seasonalized){par$ts <- 0; par$C <- 0}

    if(Ltnew < Lt.i){
      while(Lt.i > Ltnew){
        t2 <- t.i
        t1 <- t.i-tincr
        slope <- {1 - exp(-(
          par$K*(t2-t1)
          - (((par$C*par$K)/(2*pi))*sin(2*pi*(t1-par$ts)))
          + (((par$C*par$K)/(2*pi))*sin(2*pi*(t2-par$ts)))
        ))}
        Lt.i <- (par$Linf*slope - Lt.i) / (slope - 1)
        t.i <- t1
      }
      tnew <- t.i
    }

    if(Ltnew > Lt){
      while(Lt.i < Ltnew){
        t2 <- t.i+tincr
        t1 <- t.i
        slope <- {1 - exp(-(
          par$K*(t2-t1)
          - (((par$C*par$K)/(2*pi))*sin(2*pi*(t1-par$ts)))
          + (((par$C*par$K)/(2*pi))*sin(2*pi*(t2-par$ts)))
        ))}
        Lt.i <- (par$Linf*slope + Lt.i) / (slope + 1)
        t.i <- t2
      }
      tnew <- t.i
    }
  }

  # return result
  return(tnew)
}



#' @title Assign cohort, relative ages and birth date to an lfq dataset
#'
#' @description The \code{lfqCohort} function will assign cohort number (1 = oldest),
#' relative age (since length equal zero), and birth date (at length equal zero) to
#' an \code{lfq} object given a set of VBGF parameters (\code{lfq$par}).
#'
#' The method is can be used as a precursor to cohort-based analyses
#' (e.g. catch curve, VPA), or as means of identifying recruitment patterns.
#' Relative ages consider seasonalized VBGF (when applicable), thus aiding the
#' correct assignment of ages to lfq data (see Pauly 1990).
#'
#'
#' @param lfq an lfq object with fitted VBGF parameters in the \code{lfq$par} slot.
#' @param n.per.yr number of cohorts per year for slicing (Default: \code{n.per.year  = 1}).
#' See details section for more information on slicing.
#' This argument should be set to a higher values in cases of multiple cohorts per year.
#' @param agemax maximum age of the stock, which is used to define the extent of growth
#' curves via a call to \code{\link{lfqFitCurves}}. When not supplied
#' (Default: \code{agemax = NULL}), the value is estimated based on the time rquired
#' to achieve a length of 95\% of Linf.
#' @param calc_dt logical. Should time required to grow through bin be calculated.
#'   (Default: \code{calc_dt = FALSE}).
#'
#' @return a list of class `lfq` containing the additional slots for relative age
#' (\code{rel.age}), cohort number (\code{cohort}), birthday (\code{bday}), and
#' (optionally) time required to growth through bin (\code{dt}).
#'
#' @details The method involves "slicing" the lfq data based on the VBGF parameters, and
#' follows the general approach described by Pauly (1990), and demonstrated by
#' Pauly et al. (1995) using the "GOTCH.A" program. By default, a single slice per year
#' (\code{n.per.yr = 1}) is used, which sets slice boundaries of +/- 0.5 yr around the
#' \code{ta} parameter. For each lfq bin, cohort association is then determined by
#' its inclusion within a given slice. Bins of length greater than Linf, are aggregated with
#' the bin that includes Linf. In the case of multiple cohorts per year, a larger
#' \code{n.per.yr} is needed in order to resolve intra-annual cohorts.
#'
#'
#' @references
#' Pauly, D., Moreau, J., & Abad, N. (1995). Comparison of age-structured
#' and length-converted catch curves of brown trout Salmo trutta in two
#' French rivers. Fisheries Research, 22(3–4), 197–204.
#' https://doi.org/10.1016/0165-7836(94)00323-O
#'
#' Pauly, Daniel. (1990). Length-converted catch curves and the seasonal
#' growth of fishes. Fishbyte, 8(3), 33–38.
#'
#' @export
#'
#' @examples
#' # Load data and add VGBF parameters
#' data(synLFQ4)
#' synLFQ4$par <- list(Linf = 80, K = 0.5, ta = 0.25, C = 0.75, ts = 0.5)
#'
#' # visualize cohorts with lfqCohort
#' synLFQ4 <- lfqCohort(synLFQ4, n.per.yr = 1)
#' plot(synLFQ4, Fname = "catch",
#'   ylim = c(0, max(synLFQ4$midLengths)), image.col=NA)
#' pal <- colorRampPalette(c(4,5,7,2))
#' with(synLFQ4, image(x = dates, y = midLengths, z = t(cohort),
#'   col = adjustcolor(pal(max(cohort, na.rm = TRUE)), 0.75),
#'   add=TRUE
#' ))
#'
#' # visualize relative age with lfqCohort
#' plot(synLFQ4, Fname = "catch",
#'   ylim = c(0, max(synLFQ4$midLengths)), image.col=NA)
#' pal <- colorRampPalette(c(4,5,7,2))
#' with(synLFQ4, image(x = dates, y = midLengths, z = t(rel.age),
#'   col = adjustcolor(pal(100), 0.75),
#'   add=TRUE
#' ))
#'
#'
lfqCohort <- function(lfq, n.per.yr = 1,
  agemax = NULL, calc_dt = FALSE){

  if(is.null(agemax)){
    if(!is.null(lfq$agemax)){
      agemax <- lfq$agemax
    } else {
      agemax <- ceiling((1/-lfq$par$K)*log(1-((lfq$par$Linf*0.95)/lfq$par$Linf)))
    }
  }

  # record original par and make copy for adjusting ta
  PAR <- lfq$par
  PAR2 <- PAR


  ### assign cohort, rel.age, and bday for each bin
  # calc possible ta(s) (plus/minus 0.5 years to fitted lfq$par$ta)
  # defines the beginning time of the cohort ta
  tas <- sort(
    seq(
      from = (PAR$ta - 0.5),
      by = 1/n.per.yr,
      length.out = n.per.yr
    ) %% 1
  )

  # positive adjustment to agemax
  # (needed to correct for mid time of cohort)
  t_shift <- 1/n.per.yr/2

  # calc Lt for each ta; record bday to identify unique cohorts
  # this calculates the relative ages, bdays of each slice (cohort)
  Lts <- vector("list", length(tas))
  for(n in seq(Lts)){
    PAR2$ta <- tas[n]
    Lts[[n]] <- lfqFitCurves(lfq, par = PAR2,
      agemax = agemax + t_shift,
      draw = FALSE)$Lt
    Lts[[n]]$bday <- Lts[[n]]$t - Lts[[n]]$rel.age + t_shift
  }
  Lts <- do.call("rbind", Lts)

  # assign cohort number based on bday (oldest = 1)
  bdays <- sort(unique(Lts$bday))
  Lts$ct <- match(Lts$bday, bdays)

  # matching each bin to a given slice (cohort).
  # create output matrices for rel.age, cohort number, and bday
  rel.age <- cohort <- bday <- lfq$catch*NaN
  for(i in seq(length(lfq$dates))){
    t.use <- which(Lts$t == date2yeardec(lfq$dates[i]))
    Lt.use <-  Lts[t.use,]
    Lt.use <- Lt.use[order(Lt.use$Lt, decreasing = TRUE),]

    # identify upper cohort boundary
    upper <- suppressWarnings(apply(
      outer(X = lfq$midLengths, Y = Lt.use$Lt, FUN = "-"),
      MARGIN = 1,
      FUN = function(x){max(which(x < 0))}
    ))

    # If no upper exists, match to the oldest cohort (i.e. row 1 of Lt.use)
    upper <- replace(upper, upper == -Inf, 1)

    # extract result
    rel.age[,i] <- Lt.use$rel.age[upper]
    cohort[,i] <- Lt.use$ct[upper]
    bday[,i] <- Lt.use$bday[upper]
  }
  lfq$cohort <- cohort
  lfq$rel.age <- rel.age
  lfq$bday <- bday


  ### determine dt (time required to pass through bin)
  if(calc_dt){
    dt <- cohort*NaN
    dL <- diff(lfq$midLengths)/2
    lowerL <- lfq$midLengths - c(dL[1], dL)
    upperL <- lfq$midLengths + c(dL, Inf)
    for(i in seq(length(lfq$dates))){
      # t1i <- seq(length(lfq$midLengths))*NaN
      for(j in seq(length(lfq$midLengths))){
        lowert <- calc_tnew(par = PAR, Ltnew = 0,
          Lt = lowerL[j], t = lfq$dates[i], tincr = 0.01)
        uppert <- calc_tnew(par = PAR, Ltnew = 0,
          Lt = upperL[j], t = lfq$dates[i], tincr = 0.01)
        dt[j,i] <- lowert - uppert
      }
    }
    lfq$dt <- dt
  }

  return(lfq)
}



#' @title Calculate recruitment pattern of lfq object using lfqCohort approach
#'
#' @description The \code{recruitment2} function calculates the
#' recruitment pattern (i.e. monthly) of an lfq object given
#' corresponding \code{\link{VBGF}} parameters.
#' Time at recruitment is calculated via a call to
#' \code{\link{lfqCohort}}, which slices the lfq data, and assigns
#' bins to a given cohort, with associated relative ages and brithdays.
#'
#' This function does not represent a robust statistical method as it
#' assumes the same growth parameters for all individual counts
#' in the \code{lfq$catch} bins.
#' Nevertheless, the results should provide general information about
#' the recruitment pattern, e.g. relative recruitment strength by month as
#' weighted by the catch frequencies (\code{lfq$catch}).
#' Using the argument \code{use.rcounts = TRUE}, the recruitment pattern
#' is weighted by positive restructured frequencies only.
#'
#' @param lfq a length-frequency object (i.e. class \code{lfq}) that contains
#' a \code{lfq$par} slot with \code{\link{VBGF}} parameters (usually fit
#' via \code{\link{ELEFAN}}, \code{\link{ELEFAN_GA}}, or \code{\link{ELEFAN_SA}}
#' functions).
#' @param n.per.yr number of cohorts per year for slicing
#' (Default: \code{n.per.yr = 36}). Higher values may increase the resolution of the
#' recruitment estimate,  although at the expense of computation time.
#' @param agemax maximum age of the stock, which is used to define the extent of growth
#' curves via a call to \code{\link{lfqFitCurves}}. When not supplied
#' (Default: \code{agemax = NULL}), the value is estimated based on the time rquired
#' to achieve a length of 95\% of Linf.
#' @param use.rcounts logical. Should restructured counts
#' (\code{lfq$rcounts}), rather than raw catch counts (\code{lfq$catch})
#' be used to derive recruitment pattern. If true, only bins with positive restructures
#' frequencies are considered, and those values are multiplied by 100
#' in order to provide appropriate scaling for density approxmations via
#' \code{\link[graphics]{hist}}.
#' @param plot logical. Should monthly recruitment pattern [\%] by plotted
#' as a \code{\link[graphics]{barplot}} (default: \code{plot = TRUE}).
#' @param ... further arguments passed to
#' \code{\link[graphics]{barplot}} when \code{plot = TRUE}.
#'
#' @return object of "histogram" class containing
#' monthly recruitment pattern
#'
#' @details The method uses the \code{lfqCohort} function to "slice" lfq data
#' into distinct cohorts.
#' By setting the number of cohorts per year to a high value
#' (e.g. \code{n.per.yr = 36}), specific bins get associated with a precise
#' \code{ta} value, while other VBGF parameters are held constant.
#' Given that the size of recruitment can be subjective,
#' recruitment patterns should be interpreted as relative months since the
#' assumption of length zero as a recruitment size is unrealistic, as are the
#' associated precise recruitment times. Nevertheless, the overall pattern
#' (e.g. one mode vs two), should inform the resolution at which to define
#' cohorts (i.e. via \code{\link{lfqCohort}}).
#'
#' @export
#'
#' @examples
#'
#' # using catch frequencies
#' data(synLFQ4)
#' synLFQ4$par <- list(Linf = 80, K = 0.5, ta = 0.25, C = 0.75, ts = 0.5)
#' res <- recruitment(synLFQ4)
#'
#' # using restructured frequencies
#' synLFQ4 <- lfqRestructure(synLFQ4, MA = 11)
#' res <- recruitment(lfq = synLFQ4, use.rcounts = TRUE)
#'
#' # use of results in external plotting (via histogram())
#' plot(res, col = 5, freq = FALSE,
#'   main = "Recruitment pattern for synLFQ4",
#'   xlab = "Relative month"
#' )
#'
#' # presentation as percentages
#' barplot(res$counts/sum(res$counts)*100, names.arg = month.abb,
#'   las = 2, ylab = "Recruitment [%]",
#'   main = "Recruitment pattern for synLFQ4"
#' )
#'
#'
recruitment <- function(
  lfq,
  n.per.yr = 36,
  agemax = NULL,
  plot = TRUE,
  use.rcounts = FALSE,
  ...
){
  if(use.rcounts & is.null(lfq$rcounts)){stop(
    "If use.rcounts = TRUE, then lfq$rcounts can not be empty.\n
    Please run lfqRestructure first and replace lfq with the output."
  )}

  lfq <- lfqCohort(lfq = lfq, n.per.yr = n.per.yr, agemax = agemax)

  if(!use.rcounts){
    counts <- lfq$catch
  } else {
    counts <- round(lfq$rcounts*100)
    counts[which(counts<0)] <- 0
  }


  trecr <- rep(x = c(lfq$bday), times=c(counts))
  trecr <- trecr[!is.na(trecr)]
  trecr <- yeardec2date(trecr)
  trecr.mo <- as.numeric(format(trecr, "%m"))
  h <- hist(trecr.mo, breaks = seq(0.5, 12.5, by=1), plot = FALSE)
  if(plot){
    barplot(h$counts/sum(h$counts)*100,
      names.arg = 1:12,
      ylab = "Recruitment [%]",
      xlab = "Relative month",
      ...
    )
  }
  return(h)
}


#' @title Calculate age-based VBGF parameters from time-based estimates
#'
#' @description Translates ta to t0 and ts to age-adjusted ts. In addition to
#' VBGF parameters calculated via e.g. ELEFAN, the length at recruitment
#' (\code{Lrecr}) must also be provided.
#'
#' @param par list. VBGF parameters.
#' @param Lrecr numeric. Length at recruitment
#' (default: \code{Lrecr = NULL}).
#' @param plot logical. Plot graphical representation of both
#' time and age based growth curves (default: \code{plot = TRUE}).
#'
#' @return a list containging age-based VBGF parameters.
#' @export
#'
#' @examples
#' data(synLFQ4)
#' lfq <- synLFQ4
#' lfq$par <- list(
#'   Linf = 80, K = 0.5, ta = 0.25, C = 0.75, ts = 0.5
#' )
#' ta2t0(par = lfq$par, Lrecr = 10)
#'
ta2t0 <- function(par = NULL, Lrecr = NULL, plot = TRUE){
  if(is.null(par) | is.null(Lrecr)){
    stop("Error: must provide both 'par' and 'Lrecr'")
  }
  if(Lrecr > par$Linf){stop("Error: 'Lrecr' must be lower than 'par$Linf'")}

  # add seasonalized pars if needed
  seasonalized <- !is.null(par$C)
  if(!seasonalized){par$ts <- 0; par$C <- 0}

  t <- seq(0, 3, 0.01)
  Lt <- VBGF(par = par, t = t)
  trecr <- t[which.min(sqrt((Lt - Lrecr)^2))]
  t0 <- par$ta - trecr
  if(seasonalized){
    ts <- (par$ts - trecr)%%1
  } else {
    ts <- 0
  }
  par_age <- par
  par_age$ta <- NULL
  par_age$t0 <- t0
  par_age$ts <- ts
  age <- seq(t0, 3, 0.01)
  Lage <- VBGF(par = par_age, t = age)

  if(plot){
    op <- par(mfcol = c(1,2), mar = c(4,4,1,1), mgp = c(2,0.5,0), cex = 1)
    # plot by time
    plot(t, Lt, t="l", xlab = "time", ylim = c(0, max(Lt)*1.1), yaxs="i", ylab = bquote(L[t]))
    ylim <- par()$usr[3:4]
    grid(); abline(v = trecr, lty=3); box()
    lines(t, Lt)
    points(trecr, Lrecr, pch = 1, col = 2)
    points(trecr + t0, 0, pch = 4, col = 4)
    tmp <- data.frame(var = c(names(par), "trecr", "Lrecr"), val = c(unlist(par), trecr, Lrecr))
    tmp$pch = NA; tmp$col <- NA
    is_t0 <- which(tmp$var=="ta")
    tmp$pch[is_t0] <- 4
    tmp$col[is_t0] <- 4
    is_Lrecr <- which(tmp$var=="Lrecr")
    tmp$pch[is_Lrecr] <- 1
    tmp$col[is_Lrecr] <- 2
    legend("bottomright",
      legend = c(paste(tmp$var, "=", round(tmp$val,2))),
      bty = "n", cex = 0.7, pt.cex = 1,
      pch = tmp$pch, col = tmp$col
    )
    # plot by age
    plot(age, Lage, t="l", ylim = ylim, yaxs = "i", ylab = bquote(L[age]))
    grid(); abline(v = 0, lty=3); box()
    lines(age, Lage)
    points(0, Lrecr, pch = 1, col = 2)
    points(t0, 0, pch = 4, col = 4)
    tmp <- data.frame(var = c(names(par_age), "Lrecr"), val = c(unlist(par_age), Lrecr))
    tmp$pch = NA; tmp$col <- NA
    is_t0 <- which(tmp$var=="t0")
    tmp$pch[is_t0] <- 4
    tmp$col[is_t0] <- 4
    is_Lrecr <- which(tmp$var=="Lrecr")
    tmp$pch[is_Lrecr] <- 1
    tmp$col[is_Lrecr] <- 2
    legend("bottomright",
      legend = c(paste(tmp$var, "=", round(tmp$val,2))),
      bty = "n", cex = 0.7, pt.cex = 1,
      pch = tmp$pch, col = tmp$col
    )
    par(op)
  }
  return(par_age)
}

