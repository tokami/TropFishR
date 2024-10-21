#' @title Yield per recruit with selection ogive
#'
#' @description Estimates relative yield and biomass, and the first order derivative.
#'
#' @param param a list consisting of following parameters (not all are required):
#' \itemize{
#'   \item \strong{Linf}: infinite length in cm
#'   \item \strong{Winf}: infinite weight
#'   \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'   \item \strong{t0}: theoretical time zero, at which individuals of this species
#'        hatch
#'   \item \strong{M}: natural mortality
#'   \item \strong{tr}: age of recruitment
#'   \item \strong{tc}: age of first capture}
#' @param FM_change vector with ascending fishing mortalities
#' @param Lt length at time
#' @param P population size
#'
#' @keywords function prediction ypr
#'
#'
#' @details The Thompson and Bell model incorporates an iteration step
#'     simulating the stock by means of the \code{\link{stock_sim}} function. In
#'     case changes in gear characteristics - here measured in terms of Lc or
#'     tc, the length or age at first capture, respectively - should be
#'     explored, a list with selectivity information about the gear has to be
#'     provided and the prediction models make use of the selectivity
#'     \code{\link{select_ogive}} function. Sparre and Venema (1998) recommend
#'     to treat the last length class always as plus group. This model is very
#'     sensitive to zero observations in the ultimate length classes. If
#'     unrealistic results are returned, it is recommended to cut length classes
#'     with zero observations, group them in a plus group or to change the
#'     interval between length classes. Equations which are used in this
#'     function assume isometric growth, an assumption often not met. Further,
#'     the assumption that there is no relationship between the parental stock
#'     size and progeny over a wide range of fishing mortalities or exploitation
#'     values, respectively, is also said to be untrue. By default, the
#'     functions assume knife-edge recruitment and selection of gears (Sparre
#'     and Venema, 1998).
#'
#'
#' @return A list with the input parameters and dependent on the model type.
#'
#'
#' @references
#' Berkeley, S.A., and Houde, E.D., 1980. Swordfish, Xiphias gladius, dynamics in
#' the Straits of Florida. \emph{ICES C.M.}, 11.
#'
#' Beverton, R.J.H., and Holt, S.J., 1964. Table of yield functions for fishery
#' management. \emph{FAO Fish. Tech. Pap.} 38, 49 p.
#'
#' Beverton, R.J.H., and Holt, S.J., 1966. Manual of methods for fish stock
#' assessment. Pt. 2: Tables of yield functions. \emph{FAO Fisheries Technical Paper},
#' (38)Rev.1:67 p.
#'
#' Boerema, L.K., and J.A. Gulland, 1973. Stock assessment of the Peruvian anchovy
#' (Engraulis ringens) and management of the fishery. \emph{Journal of the Fisheries
#' Board of Canada}, 30(12):2226-2235
#'
#' Garcia, S. and N.P. van Zalinge, 1982. Shrimp fishing in Kuwait: methodology
#' for a joint analysis of the artisanal and industrial fisheries. pp. 119-142 In:
#' Report on the Workshop on assessment of the shrimp stocks of the west coast of
#' the Gulf between Iran and the Arabian Peninsula. Fisheries development in the
#' Gulf. Rome, FAO, FI:DP/RAB/80/015/1, 163 p.
#'
#' Gulland, J.A., 1983. Fish stock assessment: a manual of basic methods.
#' \emph{FAO/Wiley}, New York.
#'
#' Gulland, J.A. and Boerema, L., 1973. Scientific advice on catch levels.
#' \emph{Fish. Bull. (US)} 71:325-335.
#'
#' Jones, R.E. 1957. A much simplified version of the fish yield equation. Doc. No.
#' P. 21. Paper presented at the Lisbon joint meeting of International Commission
#' Northwest Atlantic-Fisheries, International Council for the Exploration of the
#' Sea, and Food and Agriculture Organization of the United Nations. 8 p. [Mimeo].
#'
#' Millar, R.B., and Holst, R., 1997. Estimation of gillnet and hook selectivity using
#' log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#' 54(3):471-477
#'
#' Pauly, D., 1980. A selection of simple methods for the assessment of tropical fish
#' stocks. \emph{FAO Fisheries Circulars (FAO)}. no. 729.
#'
#' Pauly, D., 1984. Fish population dynamics in tropical waters: a manual for use
#' with programmable calculators. \emph{ICLARM} Stud. Rev. 8, 325 p.
#'
#' Pauly, D. and M. Soriano. 1986. Some practical extensions to Beverton and
#' Holt's relative yield-per-recruit model, p. 491-495. In J.L. Maclean, L.B. Dizon
#' and L.V. Hosillos (eds.) The First Asian Fisheries Forum. Asian Fisheries Society,
#' Manila.
#'
#' Schaefer, M.B., 1954. Some aspects of the dynamics of populations important to the
#' management of the commercial marine fisheries. \emph{Inter-Am. Trop. Tuna Comm.,
#' Bull.} 1(2):27-56.
#'
#' Schaefer, M.B., 1957. A study of the dynamics of the fishery for yellowfin tuna
#' in the eastern tropical Pacific Ocean. \emph{Inter-Am. Trop. Tuna Comm., Bull.}
#' 2:247-268.
#'
#' Sparre, P., and Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.

ypr_sel <- function(param, FM_change, Lt, P){

  res <- param

  # Calculations per size class
  interval <- (Lt[2] - Lt[1])/ 2
  lower_classes <- Lt - interval
  upper_classes <- Lt + interval

  lower_t <- VBGF(param, L = lower_classes)
  upper_t <- VBGF(param, L = upper_classes)

  S1 <- (1 - (lower_classes/res$Linf))
  S2 <- (1 - (upper_classes/res$Linf))


  Z <- FM_change + res$M

  sryr <- rep(NA,length(FM_change))
  srbr <- rep(NA,length(FM_change))
  des <- rep(NA,length(FM_change))

  for(FMi in 1:length(FM_change)){
    FMx <- FM_change[FMi]
    Zx <- Z[FMi]
    # population levels
    # reduction factor per size group
    r <- (S2 ^ ((res$M/res$K) * ((FMx/Zx)/(1-(FMx/Zx)))*P)) / (S1 ^ ((res$M/res$K) * ((FMx/Zx)/(1-(FMx/Zx)))*P))

    # derivative of r
    expo <- (res$M * P * FMx) / (res$K * Zx *(1 - (FMx/Zx)))
    r.dev <- (res$M * P * S2^expo * (log(S2) - log(S1)) * Zx) / (res$K * S1^expo * (FMx-Zx)^2)

    # G per size group
    G <- rep(NA,length(Lt))
    G[1] <- r[1]  # because: rLmin-1 = 1
    for(x1 in 2:length(r)){
      G[x1] <- prod(G[x1-1], r[x1], na.rm = TRUE)
    }
    # G[length(r)] <- 0  # because: rLinf = 0

    # derivative of G
    G.dev <- rep(NA,length(Lt))
    G.dev[1] <- r.dev[1]
    for(x2a in 2:length(r)){
      G_pre <- rep(NA,x2a)
      for(x2b in 1:x2a){
        G_pre[x2b] <- r.dev[x2b] * prod(r[1:x2a][-x2b],na.rm = TRUE)
      }
      G.dev[x2a] <- sum(G_pre, na.rm = TRUE)
    }

    ry_rd1 <- ypr(param = res, FM_change = FMx, t = lower_t)
    ry_rd2 <- ypr(param = res, FM_change = FMx, t = upper_t)

    Y_R.rel_1 <- ry_rd1$ryr   #Y_R.rel_1 <- ypr.rel(FMx, Lti = lower_classes, type = "length")
    Y_R.rel_2 <- ry_rd2$ryr   # Y_R.rel_2 <- ypr.rel(FMx, Lti = upper_classes, type = "length")
    B_R.rel.tot <- ry_rd1$rbr  # B_R.rel.tot <- bpr.rel(FMx, Lti = lower_classes, type = "length")
    # derivative of Y_R.rel per size group
    dev.Y_R.rel.tot_1 <- ry_rd1$derivative  # dev.Y_R.rel.tot_1 <- derivative(FMx, Lti = lower_classes, type = "length")
    dev.Y_R.rel.tot_2 <- ry_rd2$derivative  # dev.Y_R.rel.tot_2 <- derivative(FMx, Lti = upper_classes, type = "length")

    b_1 <- rep(NA,length(S1))
    b_2 <- rep(NA,length(S1))
    for (x2 in 2:length(S1)){
      b_1[x2] <- G[x2-1] * Y_R.rel_1[x2]
      b_2[x2] <- G[x2] * Y_R.rel_2[x2]
    }
    Y_R.rel_pre <- P * (b_1-b_2)

    # derivative of b1 and b2
    dev.b_1 <- rep(NA,length(S1))
    dev.b_2 <- rep(NA,length(S1))
    for (x2 in 2:length(S1)){
      dev.b_1[x2] <- G.dev[x2-1] * dev.Y_R.rel.tot_1[x2]
      dev.b_2[x2] <- G.dev[x2] * dev.Y_R.rel.tot_2[x2]
    }
    dev.Y_R.rel_pre <- P * (dev.b_1 - dev.b_2)

    cum_Y_R.rel <- rep(0,length(S1))
    for (x3 in 2:length(S1)){
      cum_Y_R.rel[x3] <- cum_Y_R.rel[x3-1] + Y_R.rel_pre[x3]
    }

    # total yield
    nonNA <- which(!is.na(cum_Y_R.rel))
    sryr[FMi] <- cum_Y_R.rel[length(nonNA)]

    # total derivative
    dev.tot <- rep(NA,length(Lt))
    for(x5 in 2:(length(Lt)-1)){
      firstA <- Y_R.rel_pre[x5] * G.dev[x5-1] + dev.Y_R.rel_pre[x5] * G[x5-1]
      secondA <- Y_R.rel_pre[x5+1] * G.dev[x5] + dev.Y_R.rel_pre[x5+1] * G[x5]
      dev.tot[x5] <- P[x5] * (firstA - secondA)
      #dev.tot[x3] <- (P[x3]*((dev.Y_R.rel[x3]*G.dev[x3-1]) - (dev.Y_R.rel[x3+1]*G.dev[x3])))
    }
    des[FMi] <- sum(dev.tot, na.rm=TRUE)

    # total biomass
    srbr[FMi] <- B_R.rel.tot[!is.na(B_R.rel.tot)][length(B_R.rel.tot[!is.na(B_R.rel.tot)])]

  }

  ret <- c(res,list(rbr = srbr,
                    ryr = sryr,
                    derivative = des))

  return(ret)
}
