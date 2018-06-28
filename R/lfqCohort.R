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
#'
#' @return a list of class `lfq` containing the additional slots for relative age
#' (\code{rel.age}), cohort number (\code{cohort}), and birthday (\code{bday}).
#'
#' @details The method involves "slicing" the lfq data based on the VBGF parameters, and
#' follows the general approach described by Pauly (1990), and demonstrated by
#' Pauly et al. (1995) using the "GOTCH.A" program. By default, a single slice per year
#' (\code{n.per.yr = 1}) is used, which sets slice boundaries of +/- 0.5 yr around the
#' \code{t_anchor} parameter. For each lfq bin, cohort association is then determined by
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
#' synLFQ4$par <- list(Linf = 80, K = 0.5, t_anchor = 0.25, C = 0.75, ts = 0.5)
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
#' synLFQ4 <- lfqCohort(synLFQ4, n.per.yr = 1)
#' plot(synLFQ4, Fname = "catch",
#'   ylim = c(0, max(synLFQ4$midLengths)), image.col=NA)
#' pal <- colorRampPalette(c(4,5,7,2))
#' with(synLFQ4, image(x = dates, y = midLengths, z = t(rel.age),
#'   col = adjustcolor(pal(100), 0.75),
#'   add=TRUE
#' ))
#'
#'
#' ### Demonstration of catch curve on all data vs single cohort
#'
#' ### Use all data in catch curve
#' df <- data.frame(
#'   rel.age = c(synLFQ4$rel.age),
#'   n = c(synLFQ4$catch),
#'   cohort = c(synLFQ4$cohort)
#' )
#'
#' # floor of rel.age
#' df$age <- floor(df$rel.age)
#'
#' # Catch curve using age (floor rounded ages)
#' agg <- aggregate(n ~ age, data = df, FUN = sum)
#' plot(log(n) ~ age, agg)
#' fit <- lm(log(n) ~ age,
#'   data = subset(agg, age >= 1 & age <= 4)
#' )
#' abline(fit)
#' -coef(fit)[2] # true value: Z = 1.0
#'
#' # Catch curve using rel.age (unrounded ages)
#' agg <- aggregate(n ~ rel.age, data = df, FUN = sum)
#' plot(log(n) ~ rel.age, agg)
#' fit <- lm(log(n) ~ rel.age,
#'   data = subset(agg, rel.age >= 1 & rel.age <= 4)
#' )
#' abline(fit)
#' -coef(fit)[2] # true value: Z = 1.0
#'
#'
#' ### Catch curve for a single cohort
#' agg <- aggregate(n ~ cohort, df, FUN = "sum")
#' dfsub <- subset(df, cohort == 7)# most complete cohort
#'
#' # direct rel.age
#' agg <- aggregate(n ~ rel.age, data = dfsub, FUN = sum)
#' plot(log(n) ~ rel.age, agg )
#' fit <- lm(log(n) ~ rel.age,
#'   data = subset(agg, rel.age >= 1 & rel.age <= 2)
#' )
#' abline(fit)
#' -coef(fit)[2] # true value: Z = 1.0
#'
#'
#'
lfqCohort <- function(lfq, n.per.yr = 1, agemax = NULL){

  if(is.null(agemax)){
    if(!is.null(lfq$agemax)){
      agemax <- lfq$agemax
    } else {
      agemax <- ceiling((1/-lfq$par$K)*log(1-((lfq$par$Linf*0.95)/lfq$par$Linf)))
    }
  }

  # record original par and make copy for adjusting t_anchor
  PAR <- lfq$par
  PAR2 <- PAR

  # calc possible t_anchors (plus/minus 0.5 years to fitted lfq$par$t_anchor)
  t_anchors <- sort(
    seq(
      from = (PAR$t_anchor - 0.5),
      by = 1/n.per.yr,
      length.out = n.per.yr
    ) %% 1
  )

  # positive adjustment to agemax (to correct for mid time of cohort)
  t_shift <- 1/n.per.yr/2

  # calc Lt for each t_anchor; record bday to identify unique cohorts
  Lts <- vector("list", length(t_anchors))
  for(n in seq(Lts)){
    PAR2$t_anchor <- t_anchors[n]
    Lts[[n]] <- lfqFitCurves(lfq, par = PAR2,
      agemax = agemax + t_shift,
      draw = FALSE)$Lt
    Lts[[n]]$bday <- Lts[[n]]$t - Lts[[n]]$rel.age + t_shift
  }
  Lts <- do.call("rbind", Lts)

  # assign cohort number based on bday (oldest = 1)
  bdays <- sort(unique(Lts$bday))
  Lts$ct <- match(Lts$bday, bdays)

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

  return(lfq)
}
