#' @title Calculation of lengths of multiple cohorts
#'
#' @description calculates growth curves for set of growth parameters
#'
#' @param lfq a list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column),
#'   \item \strong{samplingPeriod} length of sampling period in years,
#'   \item \strong{samplingDays} time of sampling times in relation to first sampling time,
#'   \item \strong{delta_t} array with time differences between relative sampling time set to zero and
#'      other sampling times,
#'   \item \strong{rcounts} restructured frequencies,
#'   \item \strong{peaks_mat} matrix with positive peaks with distinct values,
#'   \item \strong{ASP} available sum of peaks, sum of posititve peaks which could be potential be hit by
#'      growth curves;
#' }
#' @param par a list with following growth parameters:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default: 100),
#'   \item \strong{K} curving coefficient (default: 0.1),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0.25),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; should positive peaks be flagged out?
#' @param lty The line type. Line types can either be specified as an integer
#'    (0=blank, 1=solid, 2=dashed (default), 3=dotted, 4=dotdash, 5=longdash,
#'    6=twodash) or as one of the character strings "blank", "solid", "dashed",
#'    "dotted", "dotdash", "longdash", or "twodash", where "blank" uses 'invisible
#'    lines' (i.e., does not draw them).
#' @param lwd The line width, a positive number, defaulting to 2. The interpretation
#'    is device-specific, and some devices do not implement line widths less
#'    than one. (See the help on the device for details of the interpretation.)
#' @param col A specification for the default plotting color. See section
#'    'Color Specification'.
#' @param draw logical; indicating whether growth curves should be added to
#'    existing lfq plot
#' @param tincr step for plotting
#'
#' @examples
#' data(trout)
#' res <- lfqRestructure(trout)
#' plot(res)
#' calcLt(res,par=list(Linf=16,K=0.8,t_anchor=0.28),draw=TRUE)
#'
#' @details Anchoring time point t_anchor 0.25 for example peak spawning in April or something
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{samplingPeriod}: length of sampling period in years,
#'   \item \strong{samplingDays}: time of sampling times in relation to first sampling time,
#'   \item \strong{delta_t}: array with time differences between relative sampling time set to zero and
#'      other sampling times,
#'   \item \strong{rcounts}: restructured frequencies,
#'   \item \strong{peaks_mat}: matrix with positive peaks with distinct values,
#'   \item \strong{startingPoints}: starting sample and starting length yielding in best fit;
#' }
#'
#'
#' @references
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' @export

calcLt <- function(lfq,
                   par = list(Linf = 100, K = 0.1, t_anchor = 0.25, C = 0, ts = 0),
                   agemax = NULL, flagging.out = TRUE,
                   lty = 2, lwd = 1, col = 1,
                   draw = FALSE, tincr = 0.05){

  if(is.null(par$Linf) | is.null(par$K) | is.null(par$t_anchor)) stop("Linf, K and t_anchor have to provided in par.")

  if(is.null(agemax)){
    agemax <- ceiling((1/-par$K)*log(1-((par$Linf*0.95)/par$Linf)))
  }

  yeardec <- date2yeardec(lfq$dates)

  tmax <- max(yeardec, na.rm = TRUE)
  ncohort <- agemax + ceiling(diff(range(yeardec)))
  tA.use <- par$t_anchor
  tAs <- seq(floor(tmax-ncohort), ceiling(tmax), by=1) + abs(tA.use)%%1   # tA  == tanchor
  ncohort <- length(tAs)


  par2 <- par
  t <- yeardec
  Lt <- vector(mode="list", ncohort)
  for(ct in seq(ncohort)){
    par2$t_anchor <- tAs[ct]
    Lt.ct <- VBGF(param = par2, t = t) ## do.call(what = VBGF,  par2)    # Lt.ct <- VBGF(lfq = par2, t = yeardec)   #
    Lt[[ct]] <- data.frame(t=yeardec, Lt=Lt.ct)
    if(draw){
      tmp <- par2
      t_tmp <- seq(tAs[ct], ceiling(tmax) + 1, by=tincr)
      tmp$L <- VBGF(tmp, t = t_tmp) ## do.call(what = VBGF,  tmp)
      tmp$t <- yeardec2date(t_tmp)
      lines(L ~ t, tmp, lty=lty, lwd=lwd, col=col)
    }
  }
  Lt <- do.call(rbind, Lt)


  # calc scores
  grd <- expand.grid(Lt=lfq$midLengths, t=date2yeardec(lfq$dates))
  grd$Fs <- c(lfq$rcounts)
  grd$cohort_peaks <- c(lfq$peaks_mat)
  grd$hit <- 0
  bin.width <- diff(lfq$midLengths)
  grd$bin.lower <- lfq$midLengths - c(bin.width[1], bin.width)/2
  grd$bin.upper <- lfq$midLengths + c(bin.width, bin.width[length(bin.width)])/2

  # mark all length classes (in all sampling times) which are hit
  for(h in seq(nrow(grd))){
    tmp <- grd$t[h] == Lt$t & grd$bin.lower[h] <= Lt$Lt & grd$bin.upper[h] > Lt$Lt
    if(sum(tmp, na.rm = TRUE) > 0) grd$hit[h] <-  1
    if(flagging.out){
      if(sum(tmp, na.rm = TRUE) > 0) grd$hit[h] <-  1 + grd$hit[h]
    }
  }

  # remove marks of positive peaks within one cohort peak
  if(flagging.out){
    ch <- unique(grd$cohort_peaks)
    if(0 %in% ch) ch <- ch[-which(0 %in% ch)]
    for(ci in seq(length(ch))){
      chi <- ch[ci]
      peaki <- which(grd$cohort_peaks == chi)
      dpch <- grd$hit[peaki]
      if(sum(dpch, na.rm = TRUE) > 1){
        grd$hit[peaki] <- 0
        maxi <- max(grd$Fs[peaki])
        grd$hit[((peaki[1]-1) + which(grd$Fs[peaki] == maxi))] <- 1
      }
    }
  }


  ESP <- sum(grd$Fs * grd$hit, na.rm = TRUE)
  fASP <- ESP/lfq$ASP
  fESP <- round((10^(ESP/lfq$ASP)) /10, digits = 3)
  return(list(Lt = Lt, ASP = lfq$ASP,
              ESP = ESP, fASP = fASP, fESP = fESP))
}
