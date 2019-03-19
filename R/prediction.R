#' @title Yield per recruit
#'
#' @description Estimates the absolute and relative yield and biomass per recruit and
#'    the first order derivative.
#'
#' @param pars a list consisting of following parameters (not all are required):
#'   \itemize{
#'     \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'     \item \strong{Winf}: infinite weight form a von Bertalanffy growth curve
#'    in wet weight-grams,
#'     \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'     \item \strong{ta}: time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'     \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'     \item \strong{C}: amplitude of growth oscillation of soVBGF (range: 0 to 1, default: 0),
#'     \item \strong{ts}: summer point of soVBGF (ts = WP - 0.5) (range: 0 to 1, default: 0);
#'    \item \strong{M}: natural mortality [1/year] (numeric value or vector of identical
#'      length than midLengths),
#'   \item \strong{tr}: age of recruitment,
#'   \item \strong{tc}: age of first capture,
#'   \item \strong{a}: length-weight relationship coefficent (W = a * L^b; for kg/cm3),
#'   \item \strong{b}: length-weight relationship coefficent (W = a * L^b);
#'  }
#' @param FM_change vector with ascending fishing mortalities,
#' @param t default NA
#'
#' @keywords function prediction ypr
#'
#' @details The Thompson and Bell model incorporates an iteration step simulating the
#'    stock by means
#'    of the \code{\link{stock_sim}} function. In case changes in gear
#'    characteristics -
#'    here measured in terms of Lc or tc, the length or age at first capture,
#'    respectively -
#'    should be explored, a list with selectivity information about the gear has
#'    to be provided and
#'    the prediction models make use of the selectivity \code{\link{select_ogive}}
#'    function.
#'    Sparre and Venema (1998) recommend to treat the last length class always as plus group.
#'    This model is very sensitive to zero observations in the ultimate length
#'    classes. If unrealistic results are returned,
#'    it is recommended to cut length classes with zero observations, group
#'    them in a plus group or to change the interval between
#'    length classes.
#'    Equations which are used in this function assume isometric growth, an
#'    assumption often not met. Further, the assumption that there is no relationship
#'    between the parental stock size and progeny
#'    over a wide range of fishing mortalities or exploitation values,
#'    respectively, is also said to be untrue. By default, the functions assume
#'    knife-edge recruitment and selection of gears (Sparre and Venema, 1998).
#'
#'
#' @return A list with the input parameters and dependent on the model type following
#'    list objects:
#' \itemize{
#'   \item \code{type = 'ypr'}
#' }
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

ypr <- function(pars, FM_change, t = NA){

    
    ## KNIFE EDGE
    ## t option used in sel_rel_ypr with a loop over FM

    res <- pars ## requires Winf (or Linf and a + b not implemented), K, t0, M, tr, tc

    if(is.na(t[1])) S <- exp(-res$K * (res$tc - res$t0)) ## knife edge  ## S <- (1 - (Lci_tci/Linf))
    if(!is.na(t[1])) S <- exp(-res$K * (t - res$t0)) ## for each length group for selection ogive  ##S <- (1 - (Lti/Linf))

    Z <- (FM_change + res$M)
    m <- ((1-(FM_change/Z))/(res$M/res$K))    ### == K/Z
    A <- exp(-res$M * (res$tc - res$tr))  ## A <- (((res$Linf - res$Lc)/(res$Linf - res$Lr)) ^ (res$M/res$K))

    br <- A * res$Winf * ((1/Z) - (3*S)/(Z + res$K) +
                          (3*S^2)/(Z + 2*res$K) - (S^3)/(Z + 3*res$K))

    yr <- br * FM_change

    if(is.na(t[1]))  m_p <- (1/(res$M/res$K))
    if(!is.na(t[1]))  m_p <- (1/(1-(FM_change/Z)))

    Ol <- 1-(3*S/(1+m_p))+(3*S^2/(1+(2*m_p)))-(S^3/(1+(3*m_p)))
    Ox <- (1-(FM_change/Z))*(1-((3*S)/(1+m))+((3*S^2)/(1+(2*m)))-(S^3/(1+(3*m))))

    ## rel biomass per recruit
    if(is.na(t[1])) rbr <- Ox/Ol
    if(!is.na(t[1])) rbr <- (1-(FM_change/Z)) * (Ox/Ol)

    ## rel yield per recruit
    ryr <- (FM_change/Z) * ((S)^(res$M/res$K)) * (1 - ((3*S)/(1+m)) + ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))

    ## derivative
    C <- ((res$K*(1-(FM_change/Z)))/res$M)
    B <- (S^(res$M/res$K)) * (1 - ((3*S)/(1+C)) + ((3*S^2)/(1+2*C)) - ((S^3)/(1+3*C)))
    D <- (-((3*res$K*S^3)/(res$M*((3*res$K*(1-(FM_change/Z))/res$M)+1)^2)) +
          ((6*res$K*S^2)/(res$M*((2*res$K*(1-(FM_change/Z))/res$M)+1)^2)) -
          ((3*res$K*S)/(res$M*((res$K*(1-(FM_change/Z))/res$M)+1)^2)))
    des <- (FM_change/Z) * (S^(res$M/res$K)) * D + B

    ret <- c(res, list(yr = yr,
                       br = br,
                       ryr = ryr,
                       rbr = rbr,
                       derivative = des))
    return(ret)
}

#' @title Yield per recruit with selection ogive
#'
#' @description Estimates relative yield and biomass, and the first order derivative.
#'
#' @param pars a list consisting of following parameters (not all are required):
#'   \itemize{
#'     \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'     \item \strong{Winf}: infinite weight form a von Bertalanffy growth curve
#'    in wet weight-grams.
#'     \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'     \item \strong{ta}: time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'     \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'     \item \strong{C}: amplitude of growth oscillation of soVBGF (range: 0 to 1, default: 0),
#'     \item \strong{ts}: summer point of soVBGF (ts = WP - 0.5) (range: 0 to 1, default: 0);
#'    \item \strong{M}: natural mortality [1/year] (numeric value or vector of identical
#'      length than midLengths),
#'   \item \strong{tr}: age of recruitment,
#'   \item \strong{tc}: age of first capture,
#'   \item \strong{a}: length-weight relationship coefficent (W = a * L^b; for kg/cm3),
#'   \item \strong{b}: length-weight relationship coefficent (W = a * L^b),
#'  }
#' @param FM_change vector with ascending fishing mortalities
#' @param Lt length at time
#' @param P population size
#'
#' @keywords function prediction ypr
#'
#'
#' @details The Thompson and Bell model incorporates an iteration step simulating the
#'    stock by means
#'    of the \code{\link{stock_sim}} function. In case changes in gear
#'    characteristics -
#'    here measured in terms of Lc or tc, the length or age at first capture,
#'    respectively -
#'    should be explored, a list with selectivity information about the gear has
#'    to be provided and
#'    the prediction models make use of the selectivity \code{\link{select_ogive}}
#'    function.
#'    Sparre and Venema (1998) recommend to treat the last length class always as plus group.
#'    This model is very sensitive to zero observations in the ultimate length
#'    classes. If unrealistic results are returned,
#'    it is recommended to cut length classes with zero observations, group
#'    them in a plus group or to change the interval between
#'    length classes.
#'    Equations which are used in this function assume isometric growth, an
#'    assumption often not met. Further, the assumption that there is no relationship
#'    between the parental stock size and progeny
#'    over a wide range of fishing mortalities or exploitation values,
#'    respectively, is also said to be untrue. By default, the functions assume
#'    knife-edge recruitment and selection of gears (Sparre and Venema, 1998).
#'
#'
#' @return A list with the input parameters and dependent on the model type following
#'    list objects:
#' \itemize{
#'   \item
#' }
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

ypr_sel <- function(pars, FM_change, Lt, P){

    res <- pars

    ## Calculations per size class
    interval <- (Lt[2] - Lt[1])/ 2
    lower_classes <- Lt - interval
    upper_classes <- Lt + interval

    lower_t <- VBGF(pars, L = lower_classes)
    upper_t <- VBGF(pars, L = upper_classes)

    S1 <- (1 - (lower_classes/res$Linf))
    S2 <- (1 - (upper_classes/res$Linf))


    Z <- FM_change + res$M

    sryr <- rep(NA,length(FM_change))
    srbr <- rep(NA,length(FM_change))
    des <- rep(NA,length(FM_change))

    for(FMi in 1:length(FM_change)){
        FMx <- FM_change[FMi]
        Zx <- Z[FMi]
        ## population levels
        ## reduction factor per size group
        r <- (S2 ^ ((res$M/res$K) * ((FMx/Zx)/(1-(FMx/Zx)))*P)) / (S1 ^ ((res$M/res$K) * ((FMx/Zx)/(1-(FMx/Zx)))*P))

        ## derivative of r
        expo <- (res$M * P * FMx) / (res$K * Zx *(1 - (FMx/Zx)))
        r.dev <- (res$M * P * S2^expo * (log(S2) - log(S1)) * Zx) / (res$K * S1^expo * (FMx-Zx)^2)

        ## G per size group
        G <- rep(NA,length(Lt))
        G[1] <- r[1]  ## because: rLmin-1 = 1
        for(x1 in 2:length(r)){
            G[x1] <- prod(G[x1-1], r[x1], na.rm = TRUE)
        }
        ## G[length(r)] <- 0  ## because: rLinf = 0

        ## derivative of G
        G.dev <- rep(NA,length(Lt))
        G.dev[1] <- r.dev[1]
        for(x2a in 2:length(r)){
            G_pre <- rep(NA,x2a)
            for(x2b in 1:x2a){
                G_pre[x2b] <- r.dev[x2b] * prod(r[1:x2a][-x2b],na.rm = TRUE)
            }
            G.dev[x2a] <- sum(G_pre, na.rm = TRUE)
        }

        ry_rd1 <- ypr(pars = res, FM_change = FMx, t = lower_t)
        ry_rd2 <- ypr(pars = res, FM_change = FMx, t = upper_t)

        Y_R.rel_1 <- ry_rd1$ryr   ##Y_R.rel_1 <- ypr.rel(FMx, Lti = lower_classes, type = "length")
        Y_R.rel_2 <- ry_rd2$ryr   ## Y_R.rel_2 <- ypr.rel(FMx, Lti = upper_classes, type = "length")
        B_R.rel.tot <- ry_rd1$rbr  ## B_R.rel.tot <- bpr.rel(FMx, Lti = lower_classes, type = "length")
        ## derivative of Y_R.rel per size group
        dev.Y_R.rel.tot_1 <- ry_rd1$derivative  ## dev.Y_R.rel.tot_1 <- derivative(FMx, Lti = lower_classes, type = "length")
        dev.Y_R.rel.tot_2 <- ry_rd2$derivative  ## dev.Y_R.rel.tot_2 <- derivative(FMx, Lti = upper_classes, type = "length")

        b_1 <- rep(NA,length(S1))
        b_2 <- rep(NA,length(S1))
        for (x2 in 2:length(S1)){
            b_1[x2] <- G[x2-1] * Y_R.rel_1[x2]
            b_2[x2] <- G[x2] * Y_R.rel_2[x2]
        }
        Y_R.rel_pre <- P * (b_1-b_2)

        ## derivative of b1 and b2
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

        ## total yield
        nonNA <- which(!is.na(cum_Y_R.rel))
        sryr[FMi] <- cum_Y_R.rel[length(nonNA)]

        ## total derivative
        dev.tot <- rep(NA,length(Lt))
        for(x5 in 2:(length(Lt)-1)){
            firstA <- Y_R.rel_pre[x5] * G.dev[x5-1] + dev.Y_R.rel_pre[x5] * G[x5-1]
            secondA <- Y_R.rel_pre[x5+1] * G.dev[x5] + dev.Y_R.rel_pre[x5+1] * G[x5]
            dev.tot[x5] <- P[x5] * (firstA - secondA)
            ##dev.tot[x3] <- (P[x3]*((dev.Y_R.rel[x3]*G.dev[x3-1]) - (dev.Y_R.rel[x3+1]*G.dev[x3])))
        }
        des[FMi] <- sum(dev.tot, na.rm=TRUE)

        ## total biomass
        srbr[FMi] <- B_R.rel.tot[!is.na(B_R.rel.tot)][length(B_R.rel.tot[!is.na(B_R.rel.tot)])]

    }

    ret <- c(res,list(rbr = srbr,
                      ryr = sryr,
                      derivative = des))

    return(ret)
}



#' @title Stock simulation
#' 
#' @description  This function estimates stock size, biomass and yield of a stock from
#'    fishing mortality per age class or length group. This function is embedded in the
#'    Thompson and Bell model (prediction model: \code{\link{predict_mod}}).
#'
#' @param lfq a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} or \strong{age}: midpoints of the length classes (length-frequency
#'   data) or ages (age composition data),
#'   \item \strong{meanWeight}: mean weight in kg per age class or length group,
#'   \item \strong{meanValue}: mean value per kg fish per age class or length group,
#'  \item \strong{par}: a list with growth paramters:
#'  \itemize{
#'    \item \strong{M} or \strong{Z}: natural or total mortality
#'     [1/year],
#'    \item \strong{FM}: fishing mortality [1/year] per length class
#'     (numeric value or vector of identical length than midLengths);
#'     }
#' }
#' @param age_unit indicates if the age groups are per month (\code{"month"}) or
#'    per year (\code{"year"}). Default: \code{"year"}
#' @param stock_size_1 stock size of smallest age/length group
#' @param plus_group indicates age/length group, which should be turned into a plus
#'    group (i.e. all groups above are comprised in one group)
#'
#' @keywords function prediction ypr
#'
#' @examples
#' # age-based stock simulation
#' data(shrimps)
#'
#' # option 1: without plus group
#' stock_sim(shrimps, age_unit = "month")
#'
#' # option 2: with plus group
#' stock_sim(lfq = shrimps, age_unit = "month", plus_group = 11)
#'
#' # length-based stock simulation
#' data(hake)
#'
#' stock_sim(lfq = hake, stock_size_1 = 98919.3)
#'
#' @details better to treat last group always as a plus group...
#'      if stock size 1 not provided assumes 1000 as intital population size
#'      make sure that FM is also in same unit as the classes, e.g. when classes in
#'      months than also FM has to be provided in 1/months
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{dt}: delta t,
#'   \item \strong{N}: population numbers,
#'   \item \strong{dead}: number of deaths due to natural mortality,
#'   \item \strong{C}: catch,
#'   \item \strong{Y}: yield,
#'   \item \strong{B}: biomass,
#'   \item \strong{V}: value,
#'   \item \strong{totals}: summarised output:
#'   \itemize{
#'   \item \strong{totC} total catch,
#'   \item \strong{totY} total yield,
#'   \item \strong{totV} total value,
#'   \item \strong{meanB} mean biomass.
#'   },
#' }
#'
#' @references
#' Garcia, S. and N.P. van Zalinge, 1982. Shrimp fishing in Kuwait: methodology
#'    for a joint analysis of the artisanal and industrial fisheries. pp. 119-142 In:
#'    Report on the Workshop on assessment of the shrimp stocks of the west coast of
#'    the Gulf between Iran and the Arabian Peninsula. Fisheries development in the
#'    Gulf. Rome, FAO, FI:DP/RAB/80/015/1, 163 p.
#'
#' Millar, R. B., & Holst, R. (1997). Estimation of gillnet and hook selectivity using
#' log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#' 54(3):471-477
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

stock_sim <- function(lfq, age_unit = "year",
                      stock_size_1 = NA, plus_group = NA){

    res <- lfq
    meanWeight <- res$meanWeight
    meanValue <- res$meanValue
    if(!"par" %in% names(lfq)) stop(noquote("Please provide the required parameters in res$par!"))
    par <- res$par
    

    ##mortalities
    FM <- par$FM
    if(!is.null(par$M)){
        nM <- par$M
        Z <- FM + nM
    }else{
        Z <- par$Z
        nM <- mean(Z - FM,na.rm=T)
    }

    if('midLengths' %in% names(res)) classes <- as.character(res$midLengths)
    if('age' %in% names(res)) classes <- as.character(res$age)
    ## create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    ## age based
    if('age' %in% names(res)){

        ## delta t
        dt <- c(diff(classes.num),NA)
        if(age_unit == 'month'){
            dt <- dt * 1/12
        }

        ##population per age group
        N <- rep(NA,length(classes.num))
        N[1] <- ifelse(is.na(stock_size_1),1000,stock_size_1)

        for(x1 in 2:length(N)){
            N[x1] <- N[x1-1] * exp(-Z[x1-1] * dt[x1-1])
            ## if(x1 == length(N)){
            ##   N[x1] <- N[x1-1] * exp(-Z[x1-1] * dt[x1-1])
            ## }
        }

        ##number of deaths per time step month or year
        dead <- c(abs(diff(N)),NA)

        ##catch in numbers
        C <- dead * (FM / Z)

        ##yield in kg or g
        Y <- C * meanWeight

        ##mean biomass in kg or g
        B <- Y / (FM * dt)

        ##value expressed in money units
        V <- Y * meanValue

        ##total catch, yield, value and average biomass
        totals <- data.frame(totC = sum(C, na.rm=TRUE),
                             totY = sum(Y, na.rm=TRUE),
                             totV = sum(V, na.rm=TRUE),
                             meanB = sum((B * dt), na.rm=TRUE) / sum(dt, na.rm=TRUE))   ###### more complicated biomass concept if dt is not constant, see Chapter 5

        res2 <- list(dt = dt, N = N, dead = dead, C = C,
                     Y = Y, B = B, V = V,
                     totals = totals)

        ##with plus group
        if(!is.na(plus_group)){

            df <- do.call(cbind,res2[1:7])
            df <- df[1:plus_group,]
            classes <- classes[1:plus_group]

            ##new class
            classes[length(classes)] <-
                paste(classes[length(classes)],"+",sep='')

            ##num deaths
            df[plus_group, "dead"] <- df[plus_group, "N"]

            ##catch
            new.C <- (FM[plus_group] / Z[plus_group]) * df[plus_group, "N"]
            catch.plus.dif <- new.C - df[plus_group, "C"]
            df[plus_group, "C"] <- new.C

            ##yield
            df[plus_group, "Y"] <- meanWeight[plus_group] * catch.plus.dif

            ##value
            df[plus_group, "V"] <-
                df[plus_group, "Y"] * meanValue[plus_group]

            ##biomass       ####not sure....omitted in manual
            df[plus_group, "B"] <-
                df[plus_group, "Y"] / (FM[plus_group] * df[plus_group, "dt"])


            res2 <- as.list(as.data.frame(df))

            df2 <- do.call(cbind,res)
            df2 <- df2[1:plus_group,]
            res <- as.list(as.data.frame(df2))
            res$age <- classes

            ##total catch, yield, value and average biomass
            totals <- data.frame(totC = sum(res2$C, na.rm=TRUE),
                                 totY = sum(res2$Y, na.rm=TRUE),
                                 totV = sum(res2$V, na.rm=TRUE),
                                 meanB = sum((res2$B * res2$dt), na.rm=TRUE) / sum(dt, na.rm=TRUE))   #### more complicated biomass concept if dt is not constant, see Chapter 5

            res2 <- c(res2,totals = list(totals))
        }
    }

    ## length based
    if('midLengths' %in% names(res)){

        if("par" %in% names(res)){
            Linf <- res$par$Linf
            K <- res$par$K
            ta <- ifelse("ta" %in% names(res$par), res$par$ta, 0)            
            t0 <- ifelse("t0" %in% names(res$par), res$par$t0, 0)
            C <- ifelse("C" %in% names(res$par), res$par$C, 0)
            ts <- ifelse("ts" %in% names(res$par), res$par$ts, 0)            
        }else{
            Linf <- res$Linf
            K <- res$K
            ta <- ifelse("ta" %in% names(res), res$ta, 0)            
            t0 <- ifelse("t0" %in% names(res), res$t0, 0)
            C <- ifelse("C" %in% names(res), res$C, 0)
            ts <- ifelse("ts" %in% names(res), res$ts, 0)
        }
        if(!("a" %in% names(par)) | !("b" %in% names(par))) stop("stock_sim requires information about the length-weight relationship. Please provide 'a' and 'b' estimates in res$par.")
        a <- par$a
        b <- par$b

        ##calculate size class interval
        int <- classes.num[2] - classes.num[1]

        ## t of lower and upper length classes
        lowL <- classes.num - (int / 2)
        upL <- classes.num + (int/2)

        ## H
        H <- ((Linf - lowL)/(Linf - upL)) ^ (nM/(2*K))

        ##population
        N <- rep(NA,length(classes.num))
        N[1] <- ifelse(is.na(stock_size_1), 1000, stock_size_1)
        for(x1 in 2:length(classes.num)){
            N[x1] <- N[x1-1] * ((1/H[x1-1]) - (FM[x1-1]/Z[x1-1])) /
                (H[x1-1] - (FM[x1-1]/Z[x1-1]))
        }

        ##number of deaths per time step month or year
        dead <- c(abs(diff(N)),NA)

        ##catch
        C <- dead * (FM/Z)

        ##average weight
        W <- a * ((lowL + upL)/2 ) ^ b

        ##yield
        Y <- C * W

        ##Value
        V <- Y * meanValue

        ##biomass
        B <- (dead / Z ) * W

        ##last length group
        x2 <- length(classes.num)
        C[x2] <- N[x2] * FM[x2]/Z[x2]
        W[x2] <- a * ((lowL[x2] + Linf)/2 ) ^ b
        Y[x2] <- C[x2] * W[x2]
        B[x2] <- N[x2] / Z[x2] * W[x2]
        V[x2] <- Y[x2] * meanValue[x2]

        ##total catch, yield, value and average biomass
        totals <- data.frame(totC = sum(C, na.rm=TRUE),
                             totY = sum(Y, na.rm=TRUE),
                             totV = sum(V, na.rm=TRUE),
                             meanB = sum((B), na.rm=TRUE))

        res2 <- list(N = N,
                     dead = dead,
                     C = C,
                     Y = Y,
                     B = B,
                     V = V,
                     totals = totals)
    }

    ret <- c(res,res2)
    return(ret)
}


#' @title Prediction models
#'
#' @description  This function applies Beverton & Holt's yield per recruit model
#'    as well as the Thompson & Bell model. These models predict catch, yield, biomass
#'    and economic values for different
#'    fishing mortality scenarions (in combination with gear changes).
#'
#' @param lfq a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} or \strong{age}: midpoints of the length classes (length-frequency
#'   data) or ages (age composition data; for the Thompson and Bell model),
#' \item \strong{meanWeight}: vector with mean weight per length group or age class (for the Thompson and Bell model),
#' \item \strong{meanValue}: vector with mean value per length group or age class (for the Thompson and Bell model),
#'  \item \strong{par}: a list with growth paramters:
#'  \itemize{
#'  \item \strong{Linf} or \strong{Winf}: infinite length or weight, respectively,
#'      for investigated species in cm [cm],
#'   \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'   \item \strong{t0}: theoretical time zero, at which individuals of this species
#'        hatch,
#'   \item \strong{ta}: time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'  \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'     \item \strong{C}: amplitude of growth oscillation of soVBGF (range: 0 to 1, default: 0),
#'     \item \strong{ts}: summer point of soVBGF (ts = WP - 0.5) (range: 0 to 1, default: 0);
#'     \item \strong{M}: natural mortality [1/year] (numeric value or vector of identical
#'      length than midLengths),
#'   \item \strong{Z}: total mortality,
#'   \item \strong{FM}: fishing mortality,
#'   \item \strong{a}: length-weight relationship coefficent (W = a * L^b; for kg/cm3),
#'   \item \strong{b}: length-weight relationship coefficent (W = a * L^b),
#'   \item \strong{Lr} or \strong{tr}: length or age of recruitment,
#'   \item \strong{Lc} or \strong{tc}: length or age at 50% selection;
#'  }
#' }
#' @param FM_change vector with ascending fishing mortalities (if \code{FM_relative}
#'    is set to TRUE, values can also be relative (only for Thompson and Bell model).
#'    Default are absolute values from 0 to 10). Or
#' @param E_change vector with ascending absolute exploitation rates;
#' @param FM_relative logical; indicating whether \code{FM_change} is relative or
#'    absolute. Default is FALSE (absolute fishing mortalities in \code{FM_change}).
#' @param Lc_change vector with ascending lengths at first capture (Lc), or
#' @param tc_change vector with ascending ages at first capture (tc)
#' @param type indicating which model should be applied: \code{"ypr"} for Beverton
#'    and Holt's yield per recruit model and \code{"ThompBell"} for the Thompson and Bell model
#' @param s_list list with selectivity parameters
#' @param stock_size_1 stock size of smallest size class, if NA values are calculated
#'    relative to a stock size of 1000 individuals
#' @param age_unit in which time unit the data is provided? "month" or "year"
#' @param plus_group if a value is provided, a plus group is created comprising this
#'    size class and all above
#' @param curr.Lc current Lc (length at first capture) if available
#' @param curr.E current exploitation rate if available
#' @param Lmin smallest length group where to start with selection ogive. Not required
#'    for "knife_edge" selection type
#' @param Lincr arbitrary length increment between length groups for estimation of
#'    selection ogive. The smaller the higher the resolution but the slower the model
#'    run. Not required for "knife_edge" selection type
#' @param plot logical; if TRUE results are displayed graphically
#' @param mark logical; if value of choosen points should be displayed in graph (default: TRUE)
#' @param hide.progressbar logical; should progressbar be displayed or hidden? (Default: FALSE)
#'
#' @keywords function prediction ypr
#'
#' @examples
#' #______________________________________
#' # Yiel Per Recruit (YPR) / Beverton and Holt's model
#' #______________________________________
#' # age structured data
#' # Nemipterus marginatus
#' threadfin <- list(par = list(Winf = 286, K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4))
#'
#' predict_mod(threadfin, FM_change = seq(0,6,0.1),
#'    tc_change = seq(0.2,1,0.2), type = 'ypr')  #where it is maximal  = MSY
#'
#' # Leiognathus spendens (Pauly, 1980)
#' ponyfish <- list(par = list(Winf = 64, K = 1, t0 = -0.2, M = 1.8, tr = 0.2))
#'
#' predict_mod(ponyfish, tc_change = c(0.2,0.3,1.0), type = 'ypr', plot=TRUE)
#'
#' #______________________________________
#' # length structured data
#' # Xiphias gladius (Berkeley and Houde, 1980)
#' swordfish <- list(par = list(Linf = 309, K = 0.0949, M = 0.18,
#'                   a = 0.0003, b = 3, Lr = 90))
#'
#' select.list <- list(selecType = 'trawl_ogive', L50 = 120, L75 = 132)
#' #swordfish$midLengths <- seq(60,300,5)
#'
#' output <- predict_mod(lfq = swordfish, Lc_change = c(100,118,150,180),
#'             s_list = select.list, type = 'ypr', Lmin = 90, Lincr = 8)
#' plot(output)
#'
#' data(hake)
#' hake$par$Lr <- 35
#' select.list <- list(selecType = 'trawl_ogive', L50 = 50, L75 = 54)
#' output <- predict_mod(lfq = hake, FM_change = seq(0,3,0.05),
#'                       Lc_change = seq(30,70,1), s_list = select.list,
#'                       type = 'ypr', plot = FALSE, curr.Lc = 50, curr.E = 0.73)
#' plot(output, type = "Isopleth", xaxis1 = "FM", yaxis1 = "Y_R.rel", mark = TRUE)
#'
#' output <- predict_mod(lfq = hake, E_change = seq(0,1,0.1),
#'                       Lc_change = seq(2,120,2), #s_list = select.list,
#'                       type = 'ypr', plot = FALSE)
#' plot(output, type = "Isopleth", xaxis1 = "E", yaxis1 = "B_R")
#'
#' #______________________________________
#' #      Thompson and Bell model
#' #______________________________________
#' # with age structured data
#' data(shrimps)
#'
#' output <- predict_mod(lfq = shrimps, FM_change = seq(0.1,20,0.1),
#'      type = "ThompBell", age_unit = "month", plot = TRUE)
#'
#' #______________________________________
#' # with length structured data
#' data(hake)
#' par(mar = c(5, 4, 4, 7))
#' predict_mod(lfq = hake,FM_change = seq(0.1,3,0.05),
#'      type = 'ThompBell', plot = TRUE)
#'
#' # create list with selectivity information
#' select.list <- list(selecType = 'trawl_ogive', L50 = 50, L75 = 55)
#'
#' output <- predict_mod(lfq = hake, FM_change = seq(0,2,0.1),
#'      Lc_change = seq(20,70,1),
#'      curr.E = 0.4, curr.Lc = 50,
#'      type = 'ThompBell', s_list = select.list)
#' plot(output, xaxis1 = "FM", yaxis_iso = "Lc", yaxis1 = "B_R", mark = TRUE)
#'
#'
#' @details The Thompson and Bell model incorporates an iteration step simulating the
#'    stock by means
#'    of the \code{\link{stock_sim}} function. In case changes in gear
#'    characteristics -
#'    here measured in terms of Lc or tc, the length or age at first capture,
#'    respectively -
#'    should be explored, a list with selectivity information about the gear has
#'    to be provided and
#'    the prediction models make use of the selectivity \code{\link{select_ogive}}
#'    function.
#'    Sparre and Venema (1998) recommend to treat the last length class always as plus group.
#'    This model is very sensitive to zero observations in the ultimate length
#'    classes. If unrealistic results are returned,
#'    it is recommended to cut length classes with zero observations, group
#'    them in a plus group or to change the interval between
#'    length classes.
#'    Equations which are used in this function assume isometric growth, an
#'    assumption often not met. Further, the assumption that there is no relationship
#'    between the parental stock size and progeny
#'    over a wide range of fishing mortalities or exploitation values,
#'    respectively, is also said to be untrue. By default, the functions assume
#'    knife-edge recruitment and selection of gears (Sparre and Venema, 1998).
#'    If E_change instead of FM_change is used the range is cut at E=0.9, because
#'    higher values of E correspond to unrealistic high values of fishing mortality.
#'    If no selectivity information is given (by use of s_list), knife edge selectivity
#'    with L50 equal to the first argument of Lc_change is assumed.
#'
#'
#' @return A list with the input parameters and dependent on the model type following
#'    list objects:
#' \itemize{
#'   \item \code{type = 'ypr'}
#'   \itemize{
#'      \item \strong{FM}: fishing mortalities,
#'      \item \strong{Lc} or \strong{tc}: lengths or ages at first capture,
#'      \item \strong{list_Lc_runs}: a list with dataframes for each Lc value:
#'      \itemize{
#'        \item \strong{FM_change}: fishing mortalities
#'        \item \strong{E}: expoitation rates
#'        \item \strong{Ty}: mean age in annual yield
#'        \item \strong{LY}: mean length in annual yield
#'        \item \strong{Wy}: mean weight in annual yield
#'        \item \strong{Y_R.rel}: relative yield per recruit (change in catch in
#'            weigth per recruit relative to initial Y/R value)
#'        \item \strong{B_R.rel}: relative biomass per recruit
#'        \item \strong{Y_R}: yield per recurit (catch in weight per recruit)
#'        \item \strong{B_R}: biomass per recruit
#'        \item \strong{B_R.percent}: percentage biomass per recurit in relation to virgin
#'            biomass per recruit
#'      }
#'      \item \strong{df_Es}: a dataframe with references points (columns) for
#'          different Lc values (rows)
#'      \item \strong{df_current}: a dataframe with the exploitation status, yield
#'          and biomass values of current exploitation or selectivity (if E_curr or
#'          Lc_tc_curr provided).
#'   }
#'   \item \code{type = 'ThomBell'}
#'   \itemize{
#'      \item \strong{dt}: delta t,
#'      \item \strong{N}: population number,
#'      \item \strong{dead}: deaths due to natural reasons,
#'      \item \strong{C}: catch,
#'      \item \strong{Y}: yield,
#'      \item \strong{B}: biomass,
#'      \item \strong{V}: value,
#'      \item \strong{totals}: summed up values (total catch, total yield, total value, average biomass),
#'      \item \strong{totC}: total catches for different x factors,
#'      \item \strong{totY}: total yield values for different x factors,
#'      \item \strong{totV}: total values for different x factors,
#'      \item \strong{meanB}: average biomasses for different x factors,
#'      \item \strong{F_change}: fishing mortality changes;
#'   }
#'   \item \code{type = 'ThomBell'} and \code{Lc_change} provided
#'   \itemize{
#'      \item \strong{FM_change}: fishing mortality changes,
#'      \item \strong{Lc_change}: changes in length at first capture,
#'      \item \strong{Lt}: lengths at age,
#'      \item \strong{sel}: probability of capture,
#'      \item \strong{mat_FM_Lc_com.C}: catch matrix for all fishing mortality and Lc/tc combinations,
#'      \item \strong{mat_FM_Lc_com.Y}: yield matrix for all fishing mortality and Lc/tc combinations,
#'      \item \strong{mat_FM_Lc_com.V}: value matrix for all fishing mortality and Lc/tc combinations,
#'      \item \strong{mat_FM_Lc_com.B}: biomass matrix for all fishing mortality and Lc/tc combinations;
#'   }
#' }
#'
#' @importFrom graphics plot
#' @importFrom utils setTxtProgressBar txtProgressBar
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
#'
#' @export

predict_mod <- function(lfq, type, FM_change = NA,
                        E_change = NA,
                        FM_relative = FALSE,
                        Lc_change = NULL,
                        tc_change = NULL,
                        s_list = NA,
                        stock_size_1 = NA, age_unit = 'year', curr.E = NA,
                        curr.Lc = NA,
                        plus_group = NA, Lmin = NA, Lincr = NA, plot = FALSE, mark = TRUE,
                        hide.progressbar = FALSE){
    res <- lfq
    if(!"par" %in% names(lfq)) stop(noquote("Please provide the required parameters in lfq$par!"))
    par <- res$par    
    res$FM_relative <- FM_relative

    ## Beverton and Holt's ypr
    if(type == "ypr"){
        if(FM_relative) stop(noquote("ypr does not work with relative changes in FM, please provide absolute values."))
        M <- par$M
        K <- par$K
        t0 <- ifelse(!is.null(par$t0),par$t0,0)
        a <- par$a  ## might be NULL
        b <- par$b  ## might be NULL
        Winf <- par$Winf  ## might be NULL
        Linf <- par$Linf ## might be NULL
        if("Linf" %in% names(par) & "a" %in% names(par) & "b" %in% names(par)){
            Winf <- a * (Linf ^ b)
        }
        ##Linf <- ifelse(!is.null(res$Linf),res$Linf, exp(log(Winf/a)/b))
        ## REALLY ? maybe without Linf: then message that Winf has to exist
        ##if(is.null(Linf) | is.na(Linf)) stop("Either Linf or Winf with a and b has to be provided!")
        ##if(is.null(Winf)) Winf <-  a * Linf ^ b  ######exp((log(Linf) - a)/b) # might still be NULL
        ## or              Winf <- exp(log(Linf-a)/b)

        if(length(FM_change) == 1 & is.na(FM_change[1]) &
           length(E_change) == 1 & is.na(E_change[1])){
            FM_change <- seq(0,10,0.1)
            print(noquote("No fishing mortality (FM_change) or exploitation rate (E_change) was provided, a default range for fishing mortality of 0 to 10 is used."))
        }

        ## transfer E_change into F_change if provided
        ## if(length(FM_change) == 1 & is.na(FM_change[1]) & length(E_change) != 1 & !is.na(E_change[1])){
        ##   FM_change <- (E_change * M) / (1 - E_change)
        ##   FM_change[FM_change == Inf] <- (0.9999 * M) / (1 - 0.9999)
        ## }
        if(length(FM_change) == 1 & is.na(FM_change[1]) &
           length(E_change) != 1 & !is.na(E_change[1])){
            E_change <- E_change[E_change <= 0.9]
            FM_change <- (E_change * M) / (1 - E_change)
        }

        ## Recruitment  - knife edge
        tr <- par$tr   ## might be NULL
        Lr <- par$Lr   ## might be NULL
        if(is.null(tr) & is.null(Lr)) stop("Either the age or the length at recruitment (tr or Lr) has to be provided in res$par!")
        if(!is.null(Linf)){
            if(is.null(tr)) tr <- VBGF(L=Lr,pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(L=Lr,Linf=Linf,K=K,t0=t0)
            if(is.null(Lr)) Lr <- VBGF(t=tr,pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(t=tr,Linf=Linf,K=K,t0=t0)
        }


        ## Selectivity - knife edge or with selctivtiy ogive
        tc <- par$tc   ## might be NULL
        Lc <- par$Lc   ## might be NULL
        if(is.null(tc) & is.null(Lc)){
            if("L50" %in% s_list) Lc <- s_list$L50
            if("Lc" %in% s_list) Lc <- s_list$Lc
            ##if(!("Lc" %in% s_list) & !("L50" %in% s_list))stop("Either the age or the length at first capture (tc or Lc) has to be provided in res$par! \n Or provide a Lc value in s_list!")
        }
        if(!is.null(Linf)){
            if(is.null(tc) & !is.null(Lc)) tc <- VBGF(L=Lc, pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(L=Lc,Linf=Linf,K=K,t0=t0)
            if(is.null(Lc) & !is.null(tc)) Lc <- VBGF(t=tc, pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(t=tc,Linf=Linf,K=K,t0=t0)
            if(is.null(tc_change) & !is.null(Lc_change)) tc_change <- VBGF(L=Lc_change, pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(L=Lc_change,Linf=Linf,K=K,t0=t0)
            if(is.null(Lc_change) & !is.null(tc_change)) Lc_change <- VBGF(t=tc_change, pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(t=tc_change,Linf=Linf,K=K,t0=t0)
        }
        tc <- c(tc,tc_change)
        Lc <- c(Lc,Lc_change)

        if(length(s_list) > 1){
            selecType <- s_list$selecType
        }else{
            selecType <- "knife_edge"
        }


        ## HEART
        list_Lc_runs <- vector("list", length(Lc))
        list_Es <- vector("list", length(Lc))

        ## show progress bar only if the loop has more than 1 runs
        if (!hide.progressbar) {
            nlk <- length(Lc)
            if(nlk > 1){
                pb <- txtProgressBar(min=1, max=nlk, style=3)
                counter <- 1
            }
        }

        if(is.null(Lc)) Lc_tc <- tc else Lc_tc <- Lc
        for(i in 1:length(Lc_tc)){

            Lci <- Lc[i]
            tci <- tc[i]

            Z <- (M + FM_change)
            E <- FM_change/Z

            ## KNIFE EDGE
            if(length(s_list) == 1 ){##| selecType == "knife_edge"){
                input <- list(Linf=Linf,
                              Winf = Winf,
                              K = K,
                              M = M,
                              t0 = t0,
                              tr = tr,
                              tc = tci)
                output <- ypr(input, FM_change)
                B_R <- output$br
                Y_R <- output$yr
                B_R.rel <- output$rbr
                Y_R.rel <- output$ryr
                deri <- output$derivative
            }

            ## SELECTION OGIVE
            if(length(s_list) > 1 ){##& selecType != "knife_edge"){
                if("midLengths" %in% names(res)){
                    classes <- as.character(res$midLengths)
                    ## create column without plus group (sign) if present
                    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
                    classes.num <- as.numeric(classes.num[,1])
                    Lt <- classes.num
                    ## Lincr <- res$midLengths[2] - res$midLengths[1]
                    ## Lmin <- res$midLengths[1] - (Lincr/2)
                }
                if(!"midLengths" %in% names(res)){
                    if(is.na(Lmin) | is.na(Lincr)) writeLines("No midpoints of length classes are provided. This can be done using Lmin and Lincr or by a \n midLengths element in res. A standard range from 1cm to Linf * 0.98 by 2cm is assumed")
                    Lmin <- ifelse(is.na(Lmin), 1, Lmin)
                    Lincr <- ifelse(is.na(Lincr), 2, Lincr)
                    ## mid length vector
                    Lt <- seq(Lmin,(Linf*0.98),Lincr)
                } ## make Lt fixed except depending on Linf and as detailed as possible, e.g. Lincr = 0.5 (takes a long time)

                ## selectivity
                P <- select_ogive(s_list, Lt =  Lt, Lc = Lci)

                input <- list(Linf = ifelse(length(Linf) > 0, Linf, NA),
                              Winf = ifelse(length(Winf) > 0, Winf, NA),
                              K = K,
                              M = M,
                              t0 = t0,
                              tr = tr,
                              tc = tci)
                output <- ypr_sel(input, FM_change, Lt, P)

                ## relative yield and biomass per recruit
                B_R.rel <- output$rbr
                Y_R.rel <- output$ryr
                ## derivative
                deri <- output$derivative

                ##test
                Y_R <- Y_R.rel * Winf * exp(M * (tr - t0))
                B_R <- B_R.rel * Winf * exp(M * (tr - t0))

                ## biased because only prints P for largest Lc value
                ##if(i == length(Lc)) plot(Lt, P, type = 'l', ylab = 'Prob of capture',
                ##                         main = 'Selectivity function')
            }


            ## virgin biomass
            if(0 %in% FM_change){
                Bv_R <- B_R[FM_change == 0]
            }else{
                Bv_R <- B_R[FM_change == min(FM_change,na.rm = TRUE)]
                writeLines(paste0("Biomass was not estimated for a fishing mortality (FM) of 0, thus the virgin biomass corresponds to a FM of ",min(FM_change,na.rm = TRUE)))
            }

            ##biomass in percetage of virgin biomass
            B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

            ##mean age in annual yield
            Ty <- (1 / Z) + tci

            ##mean length in the annual yield
            S <- exp(-K * (tci - t0))         ## the same:    S <- 1 - (Lci/Linf)
            Ly <- Linf * (1 - ((Z*S)/(Z+K)))

            ##mean weight in annual yield
            Wy <- (Z) * Winf *
                ((1/Z) - ((3*S)/(Z+K)) +
                 ((3*(S^2))/(Z+(2*K))) - ((S^3)/(Z + (3*K))))


            results.PBH <- data.frame(FM = FM_change,
                                      E = E)
            if(length(Ty) > 0) results.PBH$Ty <- Ty
            if(length(Ly) > 0) results.PBH$Ly <- Ly
            if(length(Wy) > 0) results.PBH$Wy <- Wy

            results.PBH$Y_R.rel <- Y_R.rel
            results.PBH$B_R.rel <- B_R.rel

            ## WHY NECESSARY???
            if(length(Y_R) > 0) results.PBH$Y_R = Y_R
            if(length(B_R) > 0) results.PBH$B_R = B_R
            if(length(B_R.percent) > 0) results.PBH$B_R.percent = B_R.percent

            list_Lc_runs[[i]] <- results.PBH


            ## reference points
            Nmsy <- which.max(Y_R.rel)  ##  should be the same as which.min(abs(deri)) which is also labelled Nmax
            deri_pot <- deri[1:Nmsy]
            N01 <- which.min(abs(deri_pot - (deri[1] * 0.1)))
            N05 <- which.min(abs(B_R.percent - 50))  ##which.min(abs(deri - (deri[1] * 0.5)))

            df_loop_Es <- data.frame(Lc = ifelse(!is.null(Lci),Lci,NA),
                                     tc = ifelse(!is.null(tci),tci,NA),
                                     F01 = FM_change[N01],
                                     Fmsy = FM_change[Nmsy])
            if(length(B_R.percent) > 0) df_loop_Es$F05 <- FM_change[N05]   ## WHY NECESSARY????
            ## df_loop_Es$Fmax <- FM_change[Nmax]
            df_loop_Es$E01 <- E[N01]
            df_loop_Es$Emsy <- E[Nmsy]
            if(length(B_R.percent) > 0) df_loop_Es$E05 <- E[N05]    ## WHY NECESSARY????
            ## df_loop_Es$Emax <- E[Nmax]

            list_Es[[i]] <- df_loop_Es

            ## update counter and progress bar
            if (!hide.progressbar) {
                if(nlk > 1){
                    setTxtProgressBar(pb, counter)
                    counter <- counter + 1
                }
            }
        }

        df_Es <- do.call(rbind,list_Es)

        names(list_Lc_runs) <- paste0("Lc_", Lc_tc)   ## names(list_tc_runs) <- tc
        ret <- c(res,list(FM_change = FM_change,
                          Lc = Lc,
                          tc = tc,
                          list_Lc_runs = list_Lc_runs,   ##   list_tc_runs = list_tc_runs,
                          df_Es = df_Es))   ##   df_Es = df_Es,


        if(!is.na(curr.E) & !is.na(curr.Lc)){
            curr.tc <- VBGF(L=curr.Lc, pars = list(Linf=Linf,K=K,t0=t0))
            ## current exploitation rate
            curr.F = (M * curr.E)/(1-curr.E)  ## curr.F <- (M * curr.E)/(1-curr.E)
            tmpList <- list(Linf=Linf,
                            Winf = Winf,
                            K = K,
                            M = M,
                            t0 = t0,
                            tr = tr,
                            tc = curr.tc)
            if(length(s_list) == 1 | selecType == "knife_edge"){
                tmpRES <- ypr(pars = tmpList, FM_change = curr.F)
            }
            if(length(s_list) > 1 & selecType != "knife_edge"){
                P <- select_ogive(s_list, Lt =  Lt, Lc = curr.Lc)
                tmpRES <- ypr_sel(pars = tmpList, FM_change = curr.F, Lt, P)
                tmpRES$yr <- tmpRES$ryr * Winf * exp(M * (tr - t0))
                tmpRES$br <- tmpRES$rbr * Winf * exp(M * (tr - t0))
            }

            df_currents <- data.frame(curr.Lc = curr.Lc,
                                      curr.tc = curr.tc,
                                      curr.E = curr.E,
                                      curr.F = curr.F,
                                      curr.YPR = tmpRES$yr,        ##ypr(curr.F, curr.Lc_tc)       ##, type = "length"),           ## curr.YPR = ypr(curr.F, curr.Lc_tc, type = "age"),
                                      curr.YPR.rel = tmpRES$ryr,     ##ypr.rel(curr.F, curr.Lc_tc),   ##, type = "length"),   ## curr.YPR.rel = ypr.rel(curr.F, curr.Lc_tc, type = "age"),
                                      curr.BPR = tmpRES$br,         ##bpr(curr.F, curr.Lc_tc),           ##, type = "length"),           ## curr.BPR = bpr(curr.F, curr.Lc_tc, type = "age"),
                                      curr.BPR.rel = tmpRES$rbr)     ##bpr.rel(curr.F, curr.Lc_tc))   ##, type = "length"))   ## curr.BPR.rel = bpr.rel(curr.F, curr.Lc_tc, type = "age"))
            ret$currents <- df_currents
        }
    }

    ## Thompson and Bell model
    if(type == "ThompBell"){
        meanWeight <- res$meanWeight
        meanValue <- res$meanValue

        ##mortalities
        FM <- par$FM
        if(is.null(par$M) & is.null(par$Z)) stop(noquote("Either M or Z (in 'res$par') has to be provided!"))
        if(!is.null(par$M)){
            nM <- par$M
            Z <- FM + nM
        }else{
            Z <- par$Z
            nM <- Z - FM
        }

        Linf <- par$Linf
        K <- par$K
        t0 <- ifelse("t0" %in% names(par),par$t0,0)

        ## Selectivity - knife edge or with selctivtiy ogive
        tc <- res$tc   ## might be NULL
        Lc <- res$Lc   ## might be NULL
        if(is.null(tc) & is.null(Lc)){
            if("L50" %in% s_list) Lc <- s_list$L50
            if("Lc" %in% s_list) Lc <- s_list$Lc
            ##if(!("Lc" %in% s_list) & !("L50" %in% s_list))stop("Either the age or the length at first capture (tc or Lc) has to be provided in res$par! \n Or provide a Lc value in s_list!")
        }
        if(!is.null(Linf)){
            if(is.null(tc) & !is.null(Lc)) tc <- VBGF(L=Lc, pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(L=Lc,Linf=Linf,K=K,t0=t0)
            if(is.null(Lc) & !is.null(tc)) Lc <- VBGF(t=tc, pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(t=tc,Linf=Linf,K=K,t0=t0)
            if(is.null(tc_change) & !is.null(Lc_change)) tc_change <- VBGF(L=Lc_change, pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(L=Lc_change,Linf=Linf,K=K,t0=t0)
            if(is.null(Lc_change) & !is.null(tc_change)) Lc_change <- VBGF(t=tc_change, pars = list(Linf=Linf,K=K,t0=t0)) ## VBGF(t=tc_change,Linf=Linf,K=K,t0=t0)
        }
        tc <- c(tc,tc_change)
        Lc <- c(Lc,Lc_change)

        ## age based
        if('age' %in% names(res)) classes <- as.character(res$age)
        ## length based
        if('midLengths' %in% names(res)) classes <- as.character(res$midLengths)

        ## create column without plus group (sign) if present
        classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
        classes.num <- as.numeric(classes.num[,1])

        if(length(FM_change) == 1 & is.na(FM_change[1]) &
           length(E_change) == 1 & is.na(E_change[1])){
            FM_change <- seq(0,10,0.1)
            print(noquote("No fishing mortality (FM_change) or exploitation rate (E_change) \nwas provided, a default range for the absolute \nfishing mortality of 0 to 10 is used."))
        }

        ## transfer E_change into F_change if provided
        if(length(FM_change) == 1 & is.na(FM_change[1]) &
           length(E_change) != 1 & !is.na(E_change[1])){
            E_change <- E_change[E_change <= 0.9]
            FM_change <- (E_change * nM) / (1 - E_change)
        }
        if(length(E_change) == 1 & is.na(E_change[1])){
            E_change <- FM_change / (FM_change + nM)
        }

        Lt <- classes.num  ## if age in names(res) Lt here is age because classes.num = age


        ## Only FM change provided without Lc_tc change
        if((is.null(tc_change) & is.null(Lc_change))){  ##  | length(s_list) == 1){

            ##if(is.null(res$par$FM) | length(res$par$FM) == 1) stop(noquote("Please provide fishing mortality FM (in 'res$par') as a vector per size class!"))

            if(is.null(res$par$FM)) stop(noquote("Please provide fishing mortality FM (in 'lfq$par')!"))
            if(length(res$par$FM) == 1){
                if(length(s_list) > 1 | !is.null(Lc[1])){
                    print(noquote("Fishing mortality per length class not povided, using selectivity information to derive fishing mortality per length class."))
                    if(length(s_list) == 1){
                        s_list <- list(selecType = "knife_edge", L50 = Lc[1])
                    }
                    sel <- select_ogive(s_list, Lt = Lt)
                    FM <- res$par$FM * sel
                }else{
                    stop(noquote("Please provide either fishing mortality FM (in 'lfq$par') per length class or a Lc value!"))
                }
            }


            ##prediction based on f_change
            if(!FM_relative){
                pred_mat <- as.matrix(FM/max(FM, na.rm = TRUE)) %*% FM_change
            }
            if(FM_relative){
                pred_mat <- as.matrix(FM) %*% FM_change
            }


            lfqX <- lfq
            pred_res_list <- list()
            for(x7 in 1:length(FM_change)){
                lfqX$par$Z <- pred_mat[,x7] + nM
                lfqX$par$FM <- pred_mat[,x7]
                resL <- stock_sim(lfqX, age_unit = age_unit,
                                  stock_size_1 = stock_size_1, plus_group = plus_group)
                pred_res_list[[x7]] <- resL$totals
            }

            pred_res_df <- do.call(rbind, pred_res_list)
            pred_res_df$FM_change <- FM_change
            pred_res_df$E_change <- E_change

            res2 <- pred_res_df
            res3 <- c(res,res2)

            ## reference points
            Bper <- rep(NA,length(pred_res_df$meanB))
            Bper[1] <- 100
            for(ix in 2:length(Bper)){
                Bper[ix] <- pred_res_df$meanB[ix]/pred_res_df$meanB[1] * 100
            }
            N05 <- which.min(abs(Bper - 50))
            Nmsy <- which.max(pred_res_df$totY)

            if(!is.null(Lc[1]) & !is.null(tc[1])){
                df_Es <- data.frame(Lc = Lc,
                                    tc = tc,
                                    Fmsy = FM_change[Nmsy],
                                    F05 = FM_change[N05],
                                    Emsy = E_change[Nmsy],
                                    E05 = E_change[N05])
            }else{
                df_Es <- data.frame(Fmsy = FM_change[Nmsy],
                                    F05 = FM_change[N05],
                                    Emsy = E_change[Nmsy],
                                    E05 = E_change[N05])
            }


            ret <- c(res3, list(df_Es = df_Es))


            if(!is.na(curr.E)){
                if(!is.na(curr.Lc)){
                    curr.tc <- VBGF(L=curr.Lc, pars = list(Linf=Linf, K=K, t0=t0))
                }else curr.tc <- NA
                ## current exploitation rate
                curr.F = (nM * curr.E)/(1-curr.E)

                if(is.na(curr.Lc)){
                    sel <- (FM / max(FM,na.rm=TRUE))
                }else if(!is.na(curr.Lc)){
                    s_list <- list(selecType = "knife_edge", L50 = curr.Lc)
                    Lt <- res$midLengths
                    sel <- select_ogive(s_list, Lt = Lt, Lc = curr.Lc)
                }
                if(length(s_list) != 1){
                    Lt <- res$midLengths
                    sel <- select_ogive(s_list, Lt = Lt)
                }

                mati <- sel * curr.F
                param.loop <- res
                param.loop$FM <- mati
                param.loop$Z <- mati + nM
                res2 <- stock_sim(lfq = param.loop, age_unit = age_unit,
                                  stock_size_1 = stock_size_1, plus_group=plus_group)
                mati2 <- res2$totals

                df_currents <- data.frame(curr.Lc = curr.Lc,
                                          curr.tc = curr.tc,
                                          curr.E = curr.E,
                                          curr.F = curr.F,
                                          curr.C = mati2$totC,
                                          curr.Y = mati2$totY,
                                          curr.V = mati2$totV,
                                          curr.B = mati2$meanB)
                ret$currents <- df_currents
            }
        }

        ## FM and Lc_tc change provided
        if(!is.null(tc_change) | !is.null(Lc_change)){
            ## instead of s_list the outcome of one of the other select functions?

            if(length(s_list) == 1){
                s_list <- list(selecType = "knife_edge", L50 = Lc[1])
            }
            sel <- select_ogive(s_list, Lt = Lt) ##classes.num

            sel.list <- list()
            for(x19 in 1:length(Lc)){
                sel.list[[x19]] <- select_ogive(s_list, Lt = Lt, Lc = Lc[x19]) ##classes.num
            }
            Lc_mat <- do.call(cbind,sel.list)
            colnames(Lc_mat) <- Lc

            Lc_mat_FM <- Lc_mat   ##max(FM, na.rm=TRUE)  ## with one it should correspond to actual fishing mortality not to change in mortality (x factor)

            ##list with FM_Lc_matrices per FM_change
            FM_Lc_com_mat.list <- list()
            if(!FM_relative){
                for(x20 in 1:length(colnames(Lc_mat_FM))){
                    FM_Lc_com_mat.list[[x20]] <- as.matrix(Lc_mat_FM[,x20]) %*% FM_change
                    colnames(FM_Lc_com_mat.list[[x20]]) <- FM_change
                }
            }
            if(FM_relative){
                for(x20 in 1:length(colnames(Lc_mat_FM))){
                    FM_Lc_com_mat.list[[x20]] <- as.matrix(Lc_mat_FM[,x20] * FM) %*% FM_change
                    colnames(FM_Lc_com_mat.list[[x20]]) <- FM_change
                }
            }

            param.loop <- res

            pred.FM_Lc_com_res_loopC_list <- vector("list",length(FM_Lc_com_mat.list))
            pred.FM_Lc_com_res_loopY_list <- vector("list",length(FM_Lc_com_mat.list))
            pred.FM_Lc_com_res_loopB_list <- vector("list",length(FM_Lc_com_mat.list))
            pred.FM_Lc_com_res_loopV_list <- vector("list",length(FM_Lc_com_mat.list))

            if (!hide.progressbar) {
                nlk <- prod(length(FM_Lc_com_mat.list),dim(FM_Lc_com_mat.list[[1]])[2])
                pb <- txtProgressBar(min=1, max=nlk, style=3)
                counter <- 1
            }

            for(x21 in 1:length(FM_Lc_com_mat.list)){  ##loop for length of list == Lc changes
                mati <- FM_Lc_com_mat.list[[x21]]

                pred.FM_Lc_com_res_loop1_list <- list()
                for(x22 in 1:dim(mati)[2]){

                    param.loop$FM <- mati[,x22]
                    param.loop$Z <- mati[,x22] + nM
                    res2 <- stock_sim(lfq = param.loop, age_unit = age_unit,
                                      stock_size_1 = stock_size_1, plus_group=plus_group)
                    pred.FM_Lc_com_res_loop1_list[[x22]] <- res2$totals

                    ## update counter and progress bar
                    if (!hide.progressbar) {
                        setTxtProgressBar(pb, counter)
                        counter <- counter + 1
                    }
                }
                prev_mat <- do.call(rbind, pred.FM_Lc_com_res_loop1_list)
                prev_matC <- prev_mat[,'totC']
                prev_matY <- prev_mat[,'totY']
                prev_matB <- prev_mat[,'meanB']
                prev_matV <- prev_mat[,'totV']

                pred.FM_Lc_com_res_loopC_list[[x21]] <- prev_matC
                pred.FM_Lc_com_res_loopY_list[[x21]] <- prev_matY
                pred.FM_Lc_com_res_loopB_list[[x21]] <- prev_matB
                pred.FM_Lc_com_res_loopV_list[[x21]] <- prev_matV
            }

            ##for catch
            mat_FM_Lc_com.C <- do.call(rbind, pred.FM_Lc_com_res_loopC_list)
            rownames(mat_FM_Lc_com.C) <- Lc
            colnames(mat_FM_Lc_com.C) <- FM_change

            ##for yield
            mat_FM_Lc_com.Y <- do.call(rbind, pred.FM_Lc_com_res_loopY_list)
            rownames(mat_FM_Lc_com.Y) <- Lc
            colnames(mat_FM_Lc_com.Y) <- FM_change

            ##for biomass
            mat_FM_Lc_com.B <- do.call(rbind, pred.FM_Lc_com_res_loopB_list)
            rownames(mat_FM_Lc_com.B) <- Lc
            colnames(mat_FM_Lc_com.B) <- FM_change

            ##for value
            mat_FM_Lc_com.V <- do.call(rbind, pred.FM_Lc_com_res_loopV_list)
            rownames(mat_FM_Lc_com.V) <- Lc
            colnames(mat_FM_Lc_com.V) <- FM_change

            ## transvers matrices for plotting (the opposite arrangement from book)
            mat_FM_Lc_com.C <- t(mat_FM_Lc_com.C)
            mat_FM_Lc_com.Y <- t(mat_FM_Lc_com.Y)
            mat_FM_Lc_com.B <- t(mat_FM_Lc_com.B)
            mat_FM_Lc_com.V <- t(mat_FM_Lc_com.V)

            ## reference points
            mat_FM_Lc_com.Bper <- matrix(NA,ncol=dim(mat_FM_Lc_com.B)[2],
                                         nrow=dim(mat_FM_Lc_com.B)[1])
            mat_FM_Lc_com.Bper[1,] <- 100
            for(ix in 2:dim(mat_FM_Lc_com.B)[1]){
                mat_FM_Lc_com.Bper[ix,] <- mat_FM_Lc_com.B[ix,]/mat_FM_Lc_com.B[1,] *100
            }
            N05 <- apply(mat_FM_Lc_com.Bper, MARGIN = 2,
                         FUN = function(x) which.min(abs(x - 50)))

            Nmsy <- apply(mat_FM_Lc_com.Y, MARGIN = 2, FUN = which.max)

            if((!is.null(Lc[1]) & !is.null(tc[1])) | (!is.na(Lc[1]) & !is.na(tc[1])) ){
                df_Es <- data.frame(Lc = Lc,
                                    tc = tc,
                                    Fmsy = FM_change[Nmsy],
                                    F05 = FM_change[N05],
                                    Emsy = E_change[Nmsy],
                                    E05 = E_change[N05])
            }else{
                df_Es <- data.frame(Fmsy = FM_change[Nmsy],
                                    F05 = FM_change[N05],
                                    Emsy = E_change[Nmsy],
                                    E05 = E_change[N05])
            }


            ret <- c(res,
                     list(FM_change = FM_change,
                          ## FM_relative = FM_relative,
                          E_change = E_change,
                          Lc_change = Lc_change,
                          tc_change = tc_change,
                          Lt = Lt,
                          sel = sel,
                          mat_FM_Lc_com.C = mat_FM_Lc_com.C,
                          mat_FM_Lc_com.Y = mat_FM_Lc_com.Y,
                          mat_FM_Lc_com.V = mat_FM_Lc_com.V,
                          mat_FM_Lc_com.B = mat_FM_Lc_com.B,
                          df_Es = df_Es))


            if(!is.na(curr.E)){
                if(!is.na(curr.Lc)){
                    curr.tc <- VBGF(L=curr.Lc, pars = list(Linf=Linf, K=K, t0=t0))
                }else curr.tc <- NA

                ## current exploitation rate
                curr.F = (nM * curr.E)/(1-curr.E)

                if(is.na(curr.Lc)){
                    sel <- FM / max(FM, na.rm = TRUE)
                }else if(!is.na(curr.Lc) | length(s_list) == 1){
                    s_list <- list(selecType = "knife_edge", L50 = curr.Lc)
                    sel <- select_ogive(s_list, Lt = Lt, Lc = curr.Lc)
                }else if(!is.na(curr.Lc) | length(s_list) != 1){
                    sel <- select_ogive(s_list, Lt = Lt, Lc = curr.Lc)
                }
                mati <- sel * curr.F
                param.loop <- res
                param.loop$FM <- mati
                param.loop$Z <- mati + nM
                res2 <- stock_sim(param.loop, age_unit,
                                  stock_size_1, plus_group=plus_group)
                mati2 <- res2$totals

                df_currents <- data.frame(curr.Lc = curr.Lc,
                                          curr.tc = curr.tc,
                                          curr.E = curr.E,
                                          curr.F = curr.F,
                                          curr.C = mati2$totC,
                                          curr.Y = mati2$totY,
                                          curr.V = mati2$totV,
                                          curr.B = mati2$meanB)
                ret$currents <- df_currents
            }
        }
    }

    
    ret$par <- par

    ## return results and plot
    class(ret) <- "predict_mod"
    if(plot) plot(ret, mark = mark)
    return(ret)
}

