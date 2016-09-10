#' @title Fitting VBGF growth curves through lfq data
#'
#' @description Thsi function estimates von Bertalanffy growth function (VBGF)
#'    curves for a set of growth parameters.
#'
#' @param lfq a list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column),
#'   \item \strong{rcounts} restructured frequencies,
#'   \item \strong{peaks_mat} matrix with positive peaks with distinct values,
#'   \item \strong{ASP} available sum of peaks, sum of posititve peaks which
#'      could be potential be hit by growth curves;
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
#'    (Default : TRUE)
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
#' data(synLFQ4)
#' res <- lfqRestructure(synLFQ4)
#' plot(res)
#' lfqFitCurves(res, par=list(Linf=80,K=0.5,t_anchor=0.75), draw=TRUE)
#'
#' @details \code{t_anchor} subsitutes the starting point from known from Fisat 2.
#'    This parameter is necessary for anchoring the growth curves on the time axis.
#'    It does not subsitute \code{t0}. However, it corresponds to the peak spawning
#'    of the species (x intercept of growth curve) and has values between 0 and 1,
#'    where 0 corresponds to spawning at the 1st of January and 0.999 corresponds to the
#'    31st of December. The default value of 0.25 or 3/12 corresponds the third month
#'    of the year, March.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{Lt}: dataframe with ages and lengths of the cohorts,
#'   \item \strong{agemax}: maximum age of species.
#'   \item \strong{ncohort}: number of cohorts,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which
#'   could be potential be hit by growth curves. This is calculated as the sum of
#'   maximum values from each run of posive restructured scores.
#'   \item \strong{ESP}: available sum of peaks,
#'   \item \strong{fASP}: available sum of peaks,
#'   \item \strong{fESP}: available sum of peaks,
#' }
#'
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a revised and expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. 1981. The relationship between gill surface area and growth performance in fish:
#' a generalization of von Bertalanffy's theory of growth. \emph{Meeresforsch}. 28:205-211
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
#' program for pocket calculators. I.C.E.S. CM 1979/6:24. Demersal Fish Cttee, 26 p.
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
#' The Collected Works of John W. Tukey-philosophy and principles of data analysis:
#' 1965-1986 (Vol. 4, pp. 517-549). Monterey, CA, USA: Wadsworth & Brooks/Cole
#'
#' @export

lfqFitCurves <- function(lfq,
                         par = list(Linf = 100, K = 0.1, t_anchor = 0.25, C = 0, ts = 0),
                         agemax = NULL, flagging.out = TRUE,
                         lty = 2, lwd = 1, col = 1,
                         draw = FALSE, tincr = 0.05){

  if(is.null(par$Linf) | is.null(par$K) | is.null(par$t_anchor)) stop("Linf, K and t_anchor have to provided in par.")

  # ISSUE: if seasonalised max age can be different and 0.95 very coarse
  if(is.null(agemax)){
    agemax <- ceiling((1/-par$K)*log(1-((par$Linf*0.95)/par$Linf)))
  }

  yeardec <- date2yeardec(lfq$dates)

  tmax <- max(yeardec, na.rm = TRUE) # maximum sample date
  ncohort <- agemax + ceiling(diff(range(yeardec))) # number of cohorts
  tA.use <- par$t_anchor # VBGF anchor time in year
  tAs <- seq(floor(tmax-ncohort), floor(tmax), by=1) + (tA.use+1e6)%%1   # anchor times
  ncohort <- length(tAs)


  par2 <- par
  t <- yeardec
  Lt <- vector(mode="list", ncohort)
  for(ct in seq(ncohort)){
    par2$t_anchor <- tAs[ct]
    rel.age <- t-tAs[ct] # relative age to anchor time
    t.ct <- t[which(rel.age <= agemax & rel.age > 0)]
    if(length(t.ct) > 0){
      Lt.ct <- VBGF(param = par2, t = t.ct) ## do.call(what = VBGF,  par2)    # Lt.ct <- VBGF(lfq = par2, t = yeardec)   #
      Lt[[ct]] <- data.frame(t=t.ct, Lt=Lt.ct)
      if(draw){
        tmp <- par2
        t_tmp <- seq(tAs[ct], max(t.ct)+tincr, by=tincr)
        tmp$L <- VBGF(tmp, t = t_tmp) ## do.call(what = VBGF,  tmp)
        tmp$t <- yeardec2date(t_tmp)
        lines(L ~ t, tmp, lty=lty, lwd=lwd, col=col)
      }
    }
  }
  Lt <- do.call(rbind, Lt)
  Lt <- Lt[which(Lt$Lt>=0),]# trim negative lengths


  # calc scores
  grd <- expand.grid(Lt=lfq$midLengths, t=date2yeardec(lfq$dates))
  grd$Fs <- c(lfq$rcounts)  # turn matrix into 1-D vector
  grd$cohort_peaks <- c(lfq$peaks_mat) # turn matrix into 1-D vector
  grd$hit <- 0 # empty vector to record bins that are "hit" by a growth curve
  bin.width <- diff(lfq$midLengths) # bin width (should allow for uneven bin sizes)
  grd$bin.lower <- lfq$midLengths - (c(bin.width[1], bin.width)/2) # upper bin limit
  grd$bin.upper <- lfq$midLengths + (c(bin.width, bin.width[length(bin.width)])/2) # lower bin limit

  # mark all length classes (in all sampling times) which are hit
  for(h in seq(nrow(grd))){
    tmp <- grd$t[h] == Lt$t & grd$bin.lower[h] <= Lt$Lt & grd$bin.upper[h] > Lt$Lt
    if(sum(tmp, na.rm = TRUE) > 0){
      if(!flagging.out) grd$hit[h] <-  1
      if(flagging.out) grd$hit[h] <- 1 + grd$hit[h]
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
        maxi <- max(grd$Fs[peaki], na.rm = TRUE)
        grd$hit[((peaki[1]-1) + which(grd$Fs[peaki] == maxi))] <- 1
      }
    }
  }


  ESP <- sum(grd$Fs * grd$hit, na.rm = TRUE)
  fASP <- ESP/lfq$ASP
  fESP <- round((10^(ESP/lfq$ASP)) /10, digits = 3)

  lfq$Lt <- Lt
  lfq$agemax <- agemax
  lfq$ncohort <- ncohort
  lfq$ESP <- ESP
  lfq$fASP <- fASP
  lfq$fESP <- fESP

  return(lfq)
}
