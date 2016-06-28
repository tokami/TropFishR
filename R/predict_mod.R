#' @title Prediction models
#'
#' @description  This function applies the Beverton & Holt's yield per recruit model
#'    as well as the Thompson & Bell model. These model predict catch, yield, biomass
#'    and economic values for different
#'    fishing mortality scenarions (in combination with gear changes).
#'
#' @param param a list consisting of following parameters (not all are required):
#' \itemize{
#'   \item \strong{Linf} or \strong{Winf}: infinite length or weight, respectively,
#'      for investigated species in cm [cm],
#'   \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'   \item \strong{t0}: theoretical time zero, at which individuals of this species
#'        hatch,
#'   \item \strong{M}: natural mortality or
#'   \item \strong{Z}: total mortality,
#'   \item \strong{FM}: fishing mortality,
#'   \item \strong{a}: length-weight relationship coefficent (W = a * L^b),
#'   \item \strong{b}: length-weight relationship coefficent (W = a * L^b),
#'   \item \strong{Lr} or \strong{tr}: length or age of recruitment;}
#' additional list objects for the Thompson and Bell model:
#'  \itemize{
#'   \item \strong{midLengths} or \strong{age}: midpoints of the length classes
#'      (length-frequency data) or ages (age composition data),
#'   \item \strong{meanWeight}: vector with mean weight per length group or age class,
#'   \item \strong{meanValue}: vector with mean value per length group or age class,
#' }
#' @param FM_change vector with ascending fishing mortalities, or
#' @param E_change vector with ascending exploitation rates
#' @param Lc_tc_change vector with ascending lengths or ages at first capture (Lc/tc)
#' @param type indicating which model should be applied: \code{"ypr"} for Beverton
#'    and Holt's yield per recruit model and \code{"ThompBell"} for the Thompson and Bell model
#' @param s_list list with selectivity parameters
#' @param stock_size_1 stock size of smallest size class, if NA values are calculated
#'    relative to a stock size of 1000 individuals
#' @param age_unit in which time unit the data is provided? "month" or "year"
#' @param plus_group if a value is provided, a plus group is created comprising this
#'    size class and all above
#' @param curr.Lc_tc current Lc (length at first capture) if available
#' @param curr.E current exploitation rate if available
#' @param Lmin smallest length group where to start with selection ogive. Not required
#'    for "knife_edge" selection type
#' @param Lincr arbitrary length increment between length groups for estimation of
#'    selection ogive. The smaller the higher the resolution but the slower the model
#'    run. Not required for "knife_edge" selection type
#' @param plot logical; if TRUE results are displayed graphically
#'
#' @keywords function prediction ypr
#'
#' @examples
#' #______________________________________
#' # Yiel Per Recruit (YPR) / Beverton and Holt's model
#' #______________________________________
#' # age structured data
#' # Nemipterus marginatus
#' threadfin <- list(Winf = 286, K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4)
#'
#' predict_mod(threadfin, FM_change = seq(0,6,0.1),
#'    Lc_tc_change = seq(0.2,1,0.2), type = 'ypr')  #where it is maximal  = MSY
#'
#' # Leiognathus spendens (Pauly, 1980)
#' ponyfish <- list(Winf = 64, K = 1, t0 = -0.2, M = 1.8, tr = 0.2)
#'
#' predict_mod(ponyfish, Lc_tc_change = c(0.2,0.3,1.0), type = 'ypr', plot=TRUE)
#'
#' #______________________________________
#' # length structured data
#' # Xiphias gladius (Berkeley and Houde, 1980)
#' swordfish <- list(Linf = 309, K = 0.0949, M = 0.18,
#'                   a = 0.0003, b = 3, Lr = 90)
#'
#' select.list <- list(selecType = 'trawl_ogive', L50 = 120, L75 = 132)
#' #swordfish$midLengths <- seq(60,300,5)
#'
#' output <- predict_mod(param = swordfish, Lc_tc_change = c(100,118,150,180),
#'             s_list = select.list, type = 'ypr', Lmin = 90, Lincr = 8)
#' plot(output)
#'
#' data(hake)
#' select.list <- list(selecType = 'trawl_ogive', L50 = 20, L75 = 24)
#' output <- predict_mod(param = hake, E_change = seq(0,1,0.05),
#'                       Lc_tc_change = seq(5,80,1), s_list = select.list,
#'                       type = 'ypr', plot = FALSE)
#' plot(output, type = "Isopleth", xaxis1 = "E", yaxis1 = "Y_R.rel", identify = FALSE)
#'
#' select.list <- list(selecType = 'knife_edge', L50 = 20, L75 = 24)
#' output <- predict_mod(param = hake, FM_change = seq(0,3,0.01),
#'                       Lc_tc_change = seq(5,80,1), s_list = select.list,
#'                       type = 'ypr', plot = FALSE)
#' plot(output, type = "Isopleth", xaxis1 = "E", yaxis1 = "B_R.rel")
#'
#' #______________________________________
#' #      Thompson and Bell model
#' #______________________________________
#' # with age structured data
#' data(shrimps)
#'
#' output <- predict_mod(param = shrimps, FM_change = seq(0.1,3,0.1),
#'      type = "ThompBell", age_unit = "month", plot = TRUE)
#'
#'
#' # create list with selectivity information
#' select.list <- list(selecType = 'trawl_ogive',
#'    L50 = 34, L75 = 36)
#'
#' # add additional parameters to data list
#' shrimps <- c(shrimps,  list(Linf = 50, K = 0.3, t0 = 0.01))
#'
#' output <- predict_mod(param = shrimps, E_change = seq(0,1,0.01),
#'    Lc_tc_change = seq(10,44,2),
#'    type = 'ThompBell', s_list = select.list,  age_unit = 'month')
#' plot(output, xaxis = "E")
#'
#' #______________________________________
#' # with length structured data
#' data(hake)
#' predict_mod(hake,FM_change = seq(0.1,3,0.1), type = 'ThompBell', plot = TRUE)
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
#'    This model is very sensitive to zeros in the ultimate length classes. If weird results are returned,
#'    it is recommended to cut length classes with zero observations or to change the interval between
#'    length classes.
#'
#'
#'    F0.1 according to Gulland and Boerema 1973
#'
#'    Caution:
#'    Equations of this function assume isometric growth (to power of 3).
#'    Often not met!
#'    Physical yield obtained by multiplying yield per recruit times number of
#'    recruits produced in the stock.
#'    Assumes that there is no relationship between parental stock size and progeny
#'    over a wide range of F or E values, respectively (which is not true). In
#'    other words: no density dependency assumed: the weight at age is constant
#'    irrespective of stock size, the natural mortality is constant irresepctive of
#'    stock size
#'    assumption of knife edge recruitment and selection implicit in equations
#'    conclusions drawn are highly dependent on the assumption of natural mortality
#'    does not consider potential population collapse at high fishing mortalities
#'    single species model
#'
#'
#' @return A list with the input parameters and dependent on the model type following
#'    list objects:
#' \itemize{
#'   \item \code{type = 'ypr'}
#'   \itemize{
#'      \item \strong{FM}: fishing mortalities,
#'      \item \strong{Lc} or \strong{tc}: lengths or ages at first capture,
#'      \item \strong{list_Lc_runs} or \strong{list_tc_runs}: a list with the dataframes for each Lc or tc value:
#'      \itemize{
#'        \item \strong{FM_change}: fishing mortality
#'        \item \strong{Y_R}: yield per recurit (catch in weight per recruit)
#'        \item \strong{Y_R.rel}: relative yield per recruit (change in catch in
#'            weigth per recruit relative to initial Y/R value)
#'        \item \strong{B_R}: biomass per recruit
#'        \item \strong{B_R.percent}: percentage biomass per recurit in relation to virgin
#'            biomass per recruit
#'        \item \strong{Ty} or \strong{LY}: mean age or mean length in annual yield
#'        \item \strong{Wy}: mean weight in annual yield
#'      }
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
#'      \item \strong{Xfact}: fishing mortality changes;
#'   }
#'   \item \code{type = 'ThomBell'} and \code{Lc_tc_change} provided
#'   \itemize{
#'      \item \strong{FM_change}: fishing mortality changes,
#'      \item \strong{Lc_tc_change}: changes in length or age at first capture,
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

predict_mod <- function(param, FM_change = NA, E_change = NA, Lc_tc_change = NULL, type,  s_list = NA,
                        stock_size_1 = NA, age_unit = 'year', curr.E = NA,
                        curr.Lc_tc = NA,
                        plus_group = NA, Lmin = NA, Lincr = NA, plot = FALSE){
  res <- param

  # Beverton and Holt's ypr
  if(type == "ypr"){
    M <- res$M
    K <- res$K
    t0 <- ifelse(!is.null(res$t0),res$t0,0)
    a <- ifelse(!is.null(res$a),res$a,NA)
    b <- ifelse(!is.null(res$b),res$b,NA)

    if(length(FM_change) == 1 & is.na(FM_change[1]) & length(E_change) == 1 & is.na(E_change[1])){
      FM_change <- seq(0,10,0.1)
      print(noquote("No fishing mortality (FM_change) or exploitation rate (E_change) was provided, a default range for fishing mortality of 0 to 10 is used."))
    }

    # transfer E_change into F_change if provided
    if(length(FM_change) == 1 & is.na(FM_change[1]) & length(E_change) != 1 & !is.na(E_change[1])){
      FM_change <- (E_change * M) / (1 - E_change)
      FM_change[FM_change == Inf] <- (0.9999 * M) / (1 - 0.9999)
    }



    # yield and biomass functions
    # type either age or length
    # ___________________________________________________
    # length yield function
    ypr <- function(FM_change, Lci_tci, type){
      if(type == "age"){
        S <- exp(-K * (Lci_tci - t0))
        A <- exp(-M * (Lci_tci - tr))
      }
      if(type == "length"){
        S <- (1 - (Lci_tci/Linf))
        A <- (((Linf - Lci_tci)/(Linf - Lr)) ^ (M/K))
      }
      Z <- (FM_change + M)
      y <- FM_change * A * Winf * ((1/Z) - (3*S)/(Z + K) +
                                     (3*S^2)/(Z + 2*K) - (S^3)/(Z + 3*K))
      return(y)
    }
    # length biomass function
    bpr <- function(FM_change, Lci_tci, type){
      if(type == "age"){
        S <- exp(-K * (Lci_tci - t0))
        A <- exp(-M * (Lci_tci - tr))
      }
      if(type == "length"){
        S <- (1 - (Lci_tci/Linf))
        A <- (((Linf - Lci_tci)/(Linf - Lr)) ^ (M/K))
      }
      Z <- (FM_change + M)
      y <- A * Winf * ((1/Z) - (3*S)/(Z + K) + (3*S^2)/(Z + 2*K) - (S^3)/(Z + 3*K))
      return(y)
    }
    # relative yield function
    ypr.rel <- function(FM_change, Lci_tci = NA, Lti = NA, type){
      if(type == "age"){
        if(!is.na(Lci_tci)) S <- exp(-K * (Lci_tci - t0)) # knife edge
        if(is.na(Lci_tci)) S <- exp(-K * (Lti - t0)) # for each length group for selection ogive
      }
      if(type == "length"){
        if(!is.na(Lci_tci)) S <- (1 - (Lci_tci/Linf))   # knife edge
        if(is.na(Lci_tci)) S <- (1 - (Lti/Linf))  # for each length group for selection ogive
      }
      Z <- (FM_change + M)
      m <- ((1-(FM_change/Z))/(M/K))    ## == K/Z
      y <- (FM_change/Z) * ((S)^(M/K)) * (1 - ((3*S)/(1+m)) + ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))
      return(y)
    }
    # relative biomass function
    bpr.rel <- function(FM_change, Lci_tci = NA, Lti = NA, type){
      if(type == "age"){
        if(!is.na(Lci_tci)) S <- exp(-K * (Lci_tci - t0)) # knife edge
        if(is.na(Lci_tci)) S <- exp(-K * (Lti - t0)) # for each length group for selection ogive
      }
      if(type == "length"){
        if(!is.na(Lci_tci)) S <- (1 - (Lci_tci/Linf))   # knife edge
        if(is.na(Lci_tci)) S <- (1 - (Lti/Linf))  # for each length group for selection ogive
      }
      Z <- (FM_change + M)
      m <- ((1-(FM_change/Z))/(M/K))    ## == K/Z
      if(!is.na(Lci_tci)) m_p <- (1/(M/K))
      if(is.na(Lci_tci)) m_p <- (1/(1-(FM_change/Z)))
      Ol <- 1-(3*S/(1+m_p))+(3*S^2/(1+(2*m_p)))-(S^3/(1+(3*m_p)))
      Ox <- (1-(FM_change/Z))*(1-((3*S)/(1+m))+((3*S^2)/(1+(2*m)))-(S^3/(1+(3*m))))
      if(!is.na(Lci_tci)) y <- Ox/Ol
      if(is.na(Lci_tci)) y <- (1-(FM_change/Z)) * (Ox/Ol)
      return(y)
    }
    # derivative of yield function
    derivative <- function(FM_change, Lci_tci = NA, Lti = NA, type){
      if(type == "age"){
        if(!is.na(Lci_tci)) S <- exp(-K * (Lci_tci - t0)) # knife edge
        if(is.na(Lci_tci)) S <- exp(-K * (Lti - t0)) # for each length group for selection ogive
      }
      if(type == "length"){
        if(!is.na(Lci_tci)) S <- (1 - (Lci_tci/Linf))   # knife edge
        if(is.na(Lci_tci)) S <- (1 - (Lti/Linf))  # for each length group for selection ogive
      }
      Z <- (FM_change + M)
      C <- ((K*(1-(FM_change/Z)))/M)
      B <- (S^(M/K)) * (1 - ((3*S)/(1+C)) + ((3*S^2)/(1+2*C)) - ((S^3)/(1+3*C)))
      D <- (-((3*K*S^3)/(M*((3*K*(1-(FM_change/Z))/M)+1)^2)) +
              ((6*K*S^2)/(M*((2*K*(1-(FM_change/Z))/M)+1)^2)) -
              ((3*K*S)/(M*((K*(1-(FM_change/Z))/M)+1)^2)))
      des <- (FM_change/Z) * (S^(M/K)) * D + B
      return(des)
    }

    # with selection ogive
    # i = length classes from Lmin to Lmax
    ypr.rel.sel <- function(FM_change, P, Lt){
      interval <- (Lt[2] - Lt[1])/ 2
      Z <- FM_change + M
      Y_R.rel.tot.all.classes <- rep(NA,length(FM_change))
      for(FMi in 1:length(FM_change)){
        FMx <- FM_change[FMi]
        Zx <- Z[FMi]
        # population levels
        # Calculations per size class
        lower_classes <- Lt - interval
        upper_classes <- Lt + interval

        S1 <- (1 - (lower_classes/Linf))
        S2 <- (1 - (upper_classes/Linf))

        # reduction factor per size group
        r <- (S2 ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P)) / (S1 ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P))

        # G per size group
        G <- rep(NA,length(Lt))
        G[1] <- r[1]  # because: rLmin-1 = 1
        for(x1 in 2:length(r)){
          G[x1] <- prod(G[x1-1], r[x1], na.rm = TRUE)
        }
        # G[length(r)] <- 0  # because: rLinf = 0

        Y_R.rel_1 <- ypr.rel(FMx, Lti = lower_classes, type = "length")
        Y_R.rel_2 <- ypr.rel(FMx, Lti = upper_classes, type = "length")

        b_1 <- rep(NA,length(S1))
        b_2 <- rep(NA,length(S1))
        for (x2 in 2:length(S1)){
          b_1[x2] <- G[x2-1] * Y_R.rel_1[x2]
          b_2[x2] <- G[x2] * Y_R.rel_2[x2]
        }

        Y_R.rel_pre <- P * (b_1-b_2)

        cum_Y_R.rel <- rep(0,length(S1))
        for (x3 in 2:length(S1)){
          cum_Y_R.rel[x3] <- cum_Y_R.rel[x3-1] + Y_R.rel_pre[x3]
        }

        nonNA <- which(!is.na(cum_Y_R.rel))
        Y_R.rel.tot.all.classes[FMi] <- cum_Y_R.rel[length(nonNA)]

      }
      return(Y_R.rel.tot.all.classes)
    }
    # relative biomass function with selction ogive
    bpr.rel.sel <- function(FM_change, P, Lt){
      Z <- FM_change + M
      B_R.rel.tot.all.classes <- rep(NA,length(FM_change))
      for(FMi in 1:length(FM_change)){

        FMx <- FM_change[FMi]
        Zx <- Z[FMi]
        # population levels
        # Calculations per size class

        S = (1 - (Lt/Linf))

        # reduction factor per size group
        r <- rep(NA, length(Lt))
        for(x1 in 2:length(Lt)){
          r[x1] <- (S[x1] ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P[x1]))  /  (S[x1-1] ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P[x1]))
        }

        # G per size group
        G <- rep(NA,length(Lt))
        for(x2 in 1:length(Lt)){
          G[x2] <- prod(r[1:x2], na.rm = TRUE)
        }
        G[1] <- r[1]  # because: rLmin-1 = 1
        G[length(Lt)] <- 0  # because: rLinf = 0


        # TWO ALTERNATIVES FOR CALCULATION NOT CLEAR:
        # 1. corresponds to Y_R / F
        B_R.rel_CLASS <- bpr.rel(FMx, Lti = Lt, type = "length")

        # 2. new approach according to Gayanilo 06
        #           m = (1-(FMx/Zx))/(M/K)  #  == K/Z
        #           mx = m / (1 - (FMx/Zx)) # == 1/(M/K)
        #           A <- (1 - ((3*S)/(1+m)) + ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))
        #           B <- (1 - ((3*S)/(1+mx)) + ((3*S^2)/(1+2*mx)) - ((S^3)/(1+3*mx)))
        #           B_R.rel_CLASS <- ((1-(FMx/Zx)) * A) / B

        B_R.rel.tot <- rep(NA,length(Lt))
        for(x3 in 2:(length(Lt)-1)){
          B_R.rel.tot[x3] <- (P[x3]*((B_R.rel_CLASS[x3]*G[x3-1]) - (B_R.rel_CLASS[x3+1]*G[x3])))
        }

        B_R.rel.tot.all.classes[FMi] <- sum(B_R.rel.tot, na.rm=TRUE)
      }
      return(B_R.rel.tot.all.classes)
    }
    bpr.rel.sel2 <- function(FM_change, P, Lt){
      interval <- (Lt[2] - Lt[1])/ 2
      Z <- FM_change + M
      B_R.rel.tot.all.classes <- rep(NA,length(FM_change))
      for(FMi in 1:length(FM_change)){

        FMx <- FM_change[FMi]
        Zx <- Z[FMi]
        # population levels
        # Calculations per size class
        lower_classes <- Lt - interval
        upper_classes <- Lt + interval

        # S1 <- (1 - (lower_classes/Linf))
        # S2 <- (1 - (upper_classes/Linf))

        # # reduction factor per size group
        # r <- (S2 ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P)) / (S1 ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P))
        #
        # # G per size group
        # G <- rep(NA,length(Lt))
        # G[1] <- r[1]  # because: rLmin-1 = 1
        # for(x1 in 2:length(r)){
        #   G[x1] <- prod(G[x1-1], r[x1], na.rm = TRUE)
        # }
        # # G[length(r)] <- 0  # because: rLinf = 0
        #
        # # TWO ALTERNATIVES FOR CALCULATION NOT CLEAR:
        # # 1. corresponds to Y_R / F
        # #B_R.rel_CLASS <- bpr.rel(FMx, Lti = Lt, type = "length")
        # B_R.rel_1 <- bpr.rel(FMx, Lti = lower_classes, type = "length")
        # B_R.rel_2 <- bpr.rel(FMx, Lti = upper_classes, type = "length")
        #
        # b_1 <- rep(NA,length(S1))
        # b_2 <- rep(NA,length(S1))
        # for (x2 in 2:length(S1)){
        #   b_1[x2] <- G[x2-1] * B_R.rel_1[x2]
        #   b_2[x2] <- G[x2] * B_R.rel_2[x2]
        # }
        #
        # B_R.rel_pre <- P * (b_1-b_2)
        #
        # cum_B_R.rel <- rep(0,length(S1))
        # for (x3 in 2:length(S1)){
        #   cum_B_R.rel[x3] <- cum_B_R.rel[x3-1] + B_R.rel_pre[x3]
        # }

        # 2. new approach according to Gayanilo 06
        #           m = (1-(FMx/Zx))/(M/K)  #  == K/Z
        #           mx = m / (1 - (FMx/Zx)) # == 1/(M/K)
        #           A <- (1 - ((3*S)/(1+m)) + ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))
        #           B <- (1 - ((3*S)/(1+mx)) + ((3*S^2)/(1+2*mx)) - ((S^3)/(1+3*mx)))
        #           B_R.rel_CLASS <- ((1-(FMx/Zx)) * A) / B

        # Sx <- (M/K)*((FMx/Zx)/(1-(FMx/Zx)))*P
        # rx <- (S2^Sx)/(S1^Sx)
        #
        # gx <- rep(NA,length(rx))
        # gx[1] <- rx[1]
        # for (x4 in 2:length(rx)) {
        #   gx[x4] <- gx[x4-1] * rx[x4]
        # }


         B_R.rel.tot <- bpr.rel(FMx, Lti = lower_classes, type = "length") ##(1-(FMx/Zx))*(N1/D1)

        #nonNA <- which(!is.na(cum_B_R.rel))
        B_R.rel.tot.all.classes[FMi] <- B_R.rel.tot[length(B_R.rel.tot)]

      }
      return(B_R.rel.tot.all.classes)
    }
    # derivative of selectivity function
    derivative.sel <- function(FM_change, P, Lt){
      Z <- (FM_change + M)
      dev.tot.all.classes <- rep(NA,length(FM_change))
      for(FMi in 1:length(FM_change)){

        FMx <- FM_change[FMi]
        Zx <- Z[FMi]
        # Calculations per size class

        S = (1 - (Lt/Linf))      ###### BIG ASSUMPTION THAT LC = Lt

        # reduction factor per size group
        # reduction factor per size group
        r <- rep(NA, length(Lt))
        for(x1 in 2:length(Lt)){
          r[x1] <- (S[x1] ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P[x1]))  /
            (S[x1-1] ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P[x1]))
        }

        # derivative of r
        r.dev <- rep(NA,length(Lt))
        for(x1a in 2:length(Lt)){
          s1 <- S[x1a]
          s2 <- S[x1a-1]
          expo <- (M * P[x1a] * FMx) / (K * Zx *(1 - (FMx/Zx)))
          r.dev[x1a] <- (M * P[x1a] * s1^expo * (log(s1) - log(s2)) * Zx) / (K * s2^expo * (FMx-Zx)^2)
        }


        # G per size group
        G <- rep(NA,length(Lt))
        for(x2 in 1:length(Lt)){
          G[x2] <- prod(r[1:x2], na.rm = TRUE)
        }
        G[1] <- r[1]  # because: rLmin-1 = 1
        G[length(Lt)] <- 0  # because: rLinf = 0

        # derivative of G
        G.dev <- rep(NA,length(Lt))
        for(x2a in 1:length(Lt)){
          re.G <- rep(NA,x2a)
          for(x2b in 1:x2a){
            if(length(r.dev[x2b] * r[1:x2a][-x2b]) == 0){
              re.G[x2b] <- NA
            }else re.G[x2b] <- r.dev[x2b] * prod(r[1:x2a][-x2b], na.rm=TRUE)
          }
          G.dev[x2a] <- sum(re.G, na.rm=TRUE)
        }

        # Yield per size group
        Y_R.rel_CLASS <- ypr.rel(FMx, Lti = Lt, type = "length")

        Y_R.rel.tot <- rep(NA,length(Lt))
        for(x3 in 2:(length(Lt)-1)){
          Y_R.rel.tot[x3] <- (P[x3]*((Y_R.rel_CLASS[x3]*G[x3-1]) - (Y_R.rel_CLASS[x3+1]*G[x3])))
        }

        # derivative of Y_R.rel per size group
        dev.Y_R.rel.tot <- derivative(FMx, Lti = Lt, type = "length")

        # total derivative
        dev.tot <- rep(NA,length(Lt))
        for(x5 in 2:(length(Lt)-1)){
          firstA <- Y_R.rel.tot[x5]*G.dev[x5-1] + dev.Y_R.rel.tot[x5]*G[x5-1]
          secondA <- Y_R.rel.tot[x5+1]*G.dev[x5] + dev.Y_R.rel.tot[x5+1]*G[x5]
          dev.tot[x5] <- P[x5] * (firstA - secondA)
          #dev.tot[x3] <- (P[x3]*((dev.Y_R.rel[x3]*G.dev[x3-1]) - (dev.Y_R.rel[x3+1]*G.dev[x3])))
        }
        dev.tot.all.classes[FMi] <- sum(dev.tot, na.rm=TRUE)
      }
      return(dev.tot.all.classes)
    }
    derivative.sel2 <- function(FM_change, P, Lt){
      interval <- (Lt[2] - Lt[1])/ 2
      Z <- (FM_change + M)
      dev.tot.all.classes <- rep(NA,length(FM_change))
      for(FMi in 1:length(FM_change)){

        FMx <- FM_change[FMi]
        Zx <- Z[FMi]
        # Calculations per size class
        lower_classes <- Lt - interval
        upper_classes <- Lt + interval

        S1 <- (1 - (lower_classes/Linf))
        S2 <- (1 - (upper_classes/Linf))

        # reduction factor per size group
        r <- (S2 ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P)) / (S1 ^ ((M/K) * ((FMx/Zx)/(1-(FMx/Zx)))*P))

        # derivative of r
        expo <- (M * P * FMx) / (K * Zx *(1 - (FMx/Zx)))
        r.dev <- (M * P * S2^expo * (log(S2) - log(S1)) * Zx) / (K * S1^expo * (FMx-Zx)^2)

        # G per size group
        G <- rep(NA,length(Lt))
        G[1] <- r[1]  # because: rLmin-1 = 1
        for(x1 in 2:length(r)){
          G[x1] <- prod(G[x1-1], r[x1], na.rm = TRUE)
        }
        # G[length(r)] <- 0  # because: rLinf = 0

        G.dev <- rep(NA,length(Lt))
        G.dev[1] <- r.dev[1]
        for(x2a in 2:length(r)){
          G_pre <- rep(NA,x2a)
          for(x2b in 1:x2a){
            G_pre[x2b] <- r.dev[x2b] * prod(r[1:x2a][-x2b],na.rm = TRUE)
          }
          G.dev[x2a] <- sum(G_pre, na.rm = TRUE)
        }

        Y_R.rel_1 <- ypr.rel(FMx, Lti = lower_classes, type = "length")
        Y_R.rel_2 <- ypr.rel(FMx, Lti = upper_classes, type = "length")

        b_1 <- rep(NA,length(S1))
        b_2 <- rep(NA,length(S1))
        for (x2 in 2:length(S1)){
          b_1[x2] <- G[x2-1] * Y_R.rel_1[x2]
          b_2[x2] <- G[x2] * Y_R.rel_2[x2]
        }

        Y_R.rel_pre <- P * (b_1-b_2)

        # derivative of Y_R.rel per size group
        dev.Y_R.rel.tot_1 <- derivative(FMx, Lti = lower_classes, type = "length")
        dev.Y_R.rel.tot_2 <- derivative(FMx, Lti = upper_classes, type = "length")

        dev.b_1 <- rep(NA,length(S1))
        dev.b_2 <- rep(NA,length(S1))
        for (x2 in 2:length(S1)){
          dev.b_1[x2] <- G.dev[x2-1] * dev.Y_R.rel.tot_1[x2]
          dev.b_2[x2] <- G.dev[x2] * dev.Y_R.rel.tot_2[x2]
        }

        dev.Y_R.rel_pre <- P * (dev.b_1 - dev.b_2)

        # total derivative
        dev.tot <- rep(NA,length(Lt))
        for(x5 in 2:(length(Lt)-1)){
          firstA <- Y_R.rel_pre[x5] * G.dev[x5-1] + dev.Y_R.rel_pre[x5] * G[x5-1]
          secondA <- Y_R.rel_pre[x5+1] * G.dev[x5] + dev.Y_R.rel_pre[x5+1] * G[x5]
          dev.tot[x5] <- P[x5] * (firstA - secondA)
          #dev.tot[x3] <- (P[x3]*((dev.Y_R.rel[x3]*G.dev[x3-1]) - (dev.Y_R.rel[x3+1]*G.dev[x3])))
        }
        dev.tot.all.classes[FMi] <- sum(dev.tot, na.rm=TRUE)

      }
      return(dev.tot.all.classes)
    }


    # ___________________________________________________

    # age based
    if("Winf" %in% names(res)){
      Winf <- res$Winf
      Linf <- a * Winf^b
      tr <- res$tr
      tc <- Lc_tc_change

      # Error message if no Lc/L50 value is provided
      if(is.na(Lc_tc_change[1]) & !("Lc" %in% s_list) & !("L50" %in% s_list)) stop("At least one Lc (L50) value has to be provided, either in Lc_tc_change or in s_list.")
      if(length(s_list) > 1){
        selecType <- s_list$selecType
      }else{
        selecType <- "knife_edge"
      }
      if(is.na(Lc_tc_change[1])){
        if("Lc" %in% s_list) sLc <- s_list$Lc
        if("L50" %in% s_list) sL50 <- s_list$L50
        Lc <- ifelse(("Lc" %in% s_list), sLc, sL50)
      }

      list_tc_runs <- vector("list",length(tc))
      list_Es <- vector("list",length(tc))

      nlk <- length(tc)
      pb <- txtProgressBar(min=1, max=nlk, style=3)
      counter <- 1

      for(i in 1:length(tc)){

        tci <- tc[i]

        Z <- (M + FM_change)
        E <- FM_change/Z

        # KNIFE EDGE
        Y_R <- ypr(FM_change, tci, type = "age")
        B_R <- bpr(FM_change, tci, type = "age")

        Y_R.rel <- ypr.rel(FM_change, tci, type = "age")
        B_R.rel <- bpr.rel(FM_change, tci, type = "age")

        deri <- derivative(FM_change, tci, type = "age")

        # SELECTION OGIVE

        #?


        #virgin biomass
        Bv_R <- B_R[which(FM_change == 0)]
        B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

        #mean age in annual yield
        Ty <- (1 / (M+FM_change)) + tci

        #mean weight in annual yield
        S <- exp(-K * (tci - t0))
        Wy <- (M+FM_change) * Winf *
          ((1/(FM_change+M)) - ((3*S)/((M+FM_change)+K)) +
             ((3*(S^2))/((M+FM_change)+(2*K))) - ((S^3)/((M+FM_change) + (3*K))))

        #mean length in the annual yield
        Ly <- Linf * (1 - (((M+FM_change)*S)/((M+FM_change)+K)))


        results.PBH <- data.frame(FM = FM_change,
                                  Ty = Ty,
                                  Wy = Wy,
                                  Ly = Ly,
                                  Y_R = Y_R,
                                  Y_R.rel = Y_R.rel,
                                  B_R = B_R,
                                  B_R.rel = B_R.rel,
                                  B_R.percent = B_R.percent)


        list_tc_runs[[i]] <- results.PBH


        # reference points
        Nmax <- which.min(abs(deri))
        deri_pot <- deri[1:Nmax]
        N01 <- which.min(abs(deri_pot - (deri[1] * 0.1)))  # this way problem that also negative part of curve can be choosen
        N05 <- which.min(abs(B_R.percent - 50)) #which.min(abs(deri_pot - (deri[1] * 0.5)))


        df_loop_Es <- data.frame(tc = tci,
                                 F01 = FM_change[N01],
                                 F05 = FM_change[N05],
                                 Fmax = FM_change[Nmax],
                                 E01 = E[N01],
                                 E05 = E[N05],
                                 Emax = E[Nmax])
        list_Es[[i]] <- df_loop_Es

        # update counter and progress bar
        setTxtProgressBar(pb, counter)
        counter <- counter + 1
      }

      df_Es <- do.call(rbind,list_Es)

      # current exploitation rate
      curr.F <- (M * curr.E)/(1-curr.E)
      df_currents <- data.frame(curr.E = curr.E,
                                curr.F = curr.F,
                                curr.YPR = ypr(curr.F, curr.Lc_tc, type = "age"),
                                curr.YPR.rel = ypr.rel(curr.F, curr.Lc_tc, type = "age"),
                                curr.BPR = bpr(curr.F, curr.Lc_tc, type = "age"),
                                curr.BPR.rel = bpr.rel(curr.F, curr.Lc_tc, type = "age"))

      names(list_tc_runs) <- tc
      ret <- c(res,list(
        FM = FM_change,
        tc = tc,
        list_tc_runs = list_tc_runs,
        df_Es = df_Es,
        currents = df_currents))
    }

    # length based
    if("Linf" %in% names(res)){

      Linf <- res$Linf
      Lr <- res$Lr
      Lc <- Lc_tc_change

      # Error message if no Lc/L50 value is provided
      if(is.na(Lc_tc_change[1]) & !("Lc" %in% s_list) & !("L50" %in% s_list)) stop("At least one Lc (L50) value has to be provided, either in Lc_tc_change or in s_list.")
      if(length(s_list) > 1){
        selecType <- s_list$selecType
      }else{
        selecType <- "knife_edge"
      }
      if(is.na(Lc_tc_change[1])){
        if("Lc" %in% s_list) sLc <- s_list$Lc
        if("L50" %in% s_list) sL50 <- s_list$L50
        Lc <- ifelse(("Lc" %in% s_list), sLc, sL50)
      }

      list_Lc_runs <- vector("list", length(Lc))
      list_Es <- vector("list", length(Lc))

      nlk <- length(Lc)
      pb <- txtProgressBar(min=1, max=nlk, style=3)
      counter <- 1

      for(i in 1:length(Lc)){

        Lci <- Lc[i]

        Z <- (M + FM_change)
        E <- FM_change/Z

        # convert Linf in Winf
        Winf <- a * (Linf ^ b)
        # convert Lr to tr
        tr <- t0 - (log(1 - (Lr / Linf)) / K )
        # convert Lc to tc
        tci <- t0 - (log(1 - (Lci / Linf)) / K )


        # KNIFE EDGE
        if(length(s_list) == 1 | selecType == "knife_edge"){
          Y_R <- ypr(FM_change, Lci, type = "length")
          B_R <- bpr(FM_change, Lci, type = "length")

          Y_R.rel <- ypr.rel(FM_change, Lci, type = "length")
          B_R.rel <- bpr.rel(FM_change, Lci, type = "length")

          deri <- derivative(FM_change, Lci, type = "length")

        }


        # SELECTION OGIVE
        if(length(s_list) > 1 & selecType != "knife_edge"){
          if("midLengths" %in% names(res)){
            classes <- as.character(res$midLengths)
            # create column without plus group (sign) if present
            classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
            classes.num <- as.numeric(classes.num[,1])
            Lt <- classes.num
            # Lincr <- res$midLengths[2] - res$midLengths[1]
            # Lmin <- res$midLengths[1] - (Lincr/2)
          }
          if(!"midLengths" %in% names(res)){
            if(is.na(Lmin) | is.na(Lincr)) stop(noquote("Please provide Lmin and Lincr!"))
            Lmin <- Lmin
            Lincr <- Lincr
            # mid length vector
            Lt <- seq(Lmin,(Linf*0.98),Lincr)
          } # make Lt fixed except depending on Linf and as detailed as possible, e.g. Lincr = 0.5 (takes a long time)

          # selectivity
          P <- select_ogive(s_list, Lt =  Lt, Lc = Lci)

          # relative yield and biomass per recruit
          Y_R.rel <- ypr.rel.sel(FM_change, P, Lt)
          B_R.rel <- bpr.rel.sel2(FM_change, P, Lt)

          #test
          Y_R <- Y_R.rel * Winf * exp(M * (tr - t0))
          B_R <- B_R.rel * Winf * exp(M * (tr - t0))

          # derivative
          deri <- derivative.sel2(FM_change, P, Lt)


          # biased because only prints P for largest Lc value
          #if(i == length(Lc)) plot(Lt, P, type = 'l', ylab = 'Prob of capture',
          #                         main = 'Selectivity function')
        }


        # virgin biomass
        if(0 %in% FM_change){
          Bv_R <- B_R[FM_change == 0]
        }else{
          Bv_R <- B_R[FM_change == min(FM_change,na.rm = TRUE)]
        }

        #biomass in percetage of virgin biomass
        B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

        #mean length in the annual yield
        S <- 1 - (Lci/Linf)
        Ly <- Linf * (1 - ((Z*S)/(Z+K)))

        #mean weight in annual yield
        Wy <- (Z) * Winf *
          ((1/Z) - ((3*S)/(Z+K)) +
             ((3*(S^2))/(Z+(2*K))) - ((S^3)/(Z + (3*K))))

        results.PBH <- data.frame(FM = FM_change,
                                  Ly = Ly,
                                  Wy = Wy,
                                  E = E,
                                  Y_R.rel = Y_R.rel,
                                  B_R.rel = B_R.rel)

        if(length(Y_R) > 0) results.PBH$Y_R = Y_R
        if(length(B_R) > 0) results.PBH$B_R = B_R
        if(length(B_R.percent) > 0) results.PBH$B_R.percent = B_R.percent


        list_Lc_runs[[i]] <- results.PBH


        # reference points
        Nmax <- which.min(abs(deri))
        deri_pot <- deri[1:Nmax]
        N01 <- which.min(abs(deri_pot - (deri[1] * 0.1)))
        N05 <- which.min(abs(B_R.percent - 50))  #which.min(abs(deri - (deri[1] * 0.5)))

        df_loop_Es <- data.frame(Lc = Lci,
                                 F01 = FM_change[N01])
        if(length(B_R.percent) > 0) df_loop_Es$F05 <- FM_change[N05]
        df_loop_Es$Fmax <- FM_change[Nmax]
        df_loop_Es$E01 <- E[N01]
        if(length(B_R.percent) > 0) df_loop_Es$E05 <- E[N05]
        df_loop_Es$Emax <- E[Nmax]

        list_Es[[i]] <- df_loop_Es

        # update counter and progress bar
        setTxtProgressBar(pb, counter)
        counter <- counter + 1
      }


      df_Es <- do.call(rbind,list_Es)

      names(list_Lc_runs) <- Lc
      ret <- c(res,list(FM = FM_change,
                        Lc = Lc,
                        list_Lc_runs = list_Lc_runs,
                        df_Es = df_Es))

      if(!is.na(curr.E) & !is.na(curr.Lc_tc)){
        # current exploitation rate
        curr.F = (M * curr.E)/(1-curr.E)
        df_currents <- data.frame(curr.E = curr.E,
                                  curr.F = curr.F,
                                  curr.YPR = ypr(curr.F, curr.Lc_tc, type = "length"),
                                  curr.YPR.rel = ypr.rel(curr.F, curr.Lc_tc, type = "length"),
                                  curr.BPR = bpr(curr.F, curr.Lc_tc, type = "length"),
                                  curr.BPR.rel = bpr.rel(curr.F, curr.Lc_tc, type = "length"))
        ret$currents <- df_currents
      }
    }
  }

  # Thompson and Bell model
  if(type == "ThompBell"){
    meanWeight <- res$meanWeight
    meanValue <- res$meanValue

    #mortalities
    FM <- res$FM
    if(is.null(res$M) & is.null(res$Z)) stop(noquote("Either M or Z (in 'param') has to be provided!"))
    if(!is.null(res$M)){
      nM <- res$M
      Z <- FM + nM
    }else{
      Z <- res$Z
      nM <- Z - FM
    }

    Linf <- res$Linf
    K <- res$K
    t0 <- ifelse("t0" %in% names(res),res$t0,0)

    # age based
    if('age' %in% names(res)) classes <- as.character(res$age)
    # length based
    if('midLengths' %in% names(res)) classes <- as.character(res$midLengths)

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    # transfer E_change into F_change if provided
    if(length(FM_change) == 1 & is.na(FM_change[1]) & length(E_change) != 1 & !is.na(E_change[1])){
      FM_change <- (E_change * nM) / (1 - E_change)
      FM_change[FM_change == Inf] <- (0.9999 * nM) / (1 - 0.9999)
    }

    # Only FM change provided without Lc_tc change
    if(is.null(Lc_tc_change) | length(s_list) == 1){
      #prediction based on f_change
      pred_mat <- as.matrix(FM) %*% FM_change

      pred_res_list <- list()
      for(x7 in 1:length(FM_change)){
        param$Z <- pred_mat[,x7] + nM
        param$FM <- pred_mat[,x7]
        res <- stock_sim(param)
        pred_res_list[[x7]] <- res$totals
      }

      pred_res_df <- do.call(rbind, pred_res_list)
      pred_res_df$Xfact <- FM_change

      res2 <- pred_res_df
      ret <- c(res,res2)
    }

    # FM and Lc_tc change provided
    if(!is.null(Lc_tc_change) & length(s_list) > 1){
      # instead of s_list the outcome of one of the other select functions?


      Lt <- Linf * (1- exp(-K * (classes.num - t0)))

      sel <- select_ogive(s_list,Lt = Lt) #classes.num

      sel.list <- list()
      for(x19 in 1:length(Lc_tc_change)){
        sel.list[[x19]] <- select_ogive(s_list, Lt = Lt, Lc = Lc_tc_change[x19]) #classes.num
      }
      Lc_mat <- do.call(cbind,sel.list)
      colnames(Lc_mat) <- Lc_tc_change

      Lc_mat_FM <- Lc_mat * max(FM, na.rm=TRUE)

      #list with FM_Lc_matrices per FM_change
      FM_Lc_com_mat.list <- list()
      for(x20 in 1:length(colnames(Lc_mat_FM))){
        FM_Lc_com_mat.list[[x20]] <- as.matrix(Lc_mat_FM[,x20]) %*% FM_change
        colnames(FM_Lc_com_mat.list[[x20]]) <- FM_change
      }

      param.loop <- res

      pred.FM_Lc_com_res_loopC_list <- vector("list",length(FM_Lc_com_mat.list))
      pred.FM_Lc_com_res_loopY_list <- vector("list",length(FM_Lc_com_mat.list))
      pred.FM_Lc_com_res_loopB_list <- vector("list",length(FM_Lc_com_mat.list))
      pred.FM_Lc_com_res_loopV_list <- vector("list",length(FM_Lc_com_mat.list))

      nlk <- prod(length(FM_Lc_com_mat.list),dim(FM_Lc_com_mat.list[[1]])[2])
      pb <- txtProgressBar(min=1, max=nlk, style=3)
      counter <- 1

      for(x21 in 1:length(FM_Lc_com_mat.list)){  #loop for length of list == Lc changes
        mati <- FM_Lc_com_mat.list[[x21]]

        pred.FM_Lc_com_res_loop1_list <- list()
        for(x22 in 1:dim(mati)[2]){

          param.loop$FM <- mati[,x22]
          param.loop$Z <- mati[,x22] + nM
          res2 <- stock_sim(param.loop, age_unit,
                            stock_size_1, plus_group=plus_group)
          pred.FM_Lc_com_res_loop1_list[[x22]] <- res2$totals

          # update counter and progress bar
          setTxtProgressBar(pb, counter)
          counter <- counter + 1
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

      #for catch
      mat_FM_Lc_com.C <- do.call(rbind, pred.FM_Lc_com_res_loopC_list)
      rownames(mat_FM_Lc_com.C) <- Lc_tc_change
      colnames(mat_FM_Lc_com.C) <- FM_change

      #for yield
      mat_FM_Lc_com.Y <- do.call(rbind, pred.FM_Lc_com_res_loopY_list)
      rownames(mat_FM_Lc_com.Y) <- Lc_tc_change
      colnames(mat_FM_Lc_com.Y) <- FM_change

      #for biomass
      mat_FM_Lc_com.B <- do.call(rbind, pred.FM_Lc_com_res_loopB_list)
      rownames(mat_FM_Lc_com.B) <- Lc_tc_change
      colnames(mat_FM_Lc_com.B) <- FM_change

      #for value
      mat_FM_Lc_com.V <- do.call(rbind, pred.FM_Lc_com_res_loopV_list)
      rownames(mat_FM_Lc_com.V) <- Lc_tc_change
      colnames(mat_FM_Lc_com.V) <- FM_change

      # transvers matrices for plotting (the opposite arrangement from book)
      mat_FM_Lc_com.C <- t(mat_FM_Lc_com.C)
      mat_FM_Lc_com.Y <- t(mat_FM_Lc_com.Y)
      mat_FM_Lc_com.B <- t(mat_FM_Lc_com.B)
      mat_FM_Lc_com.V <- t(mat_FM_Lc_com.V)

      ret <- c(res,
               list(FM_change = FM_change,
                    Lc_tc_change = Lc_tc_change,
                    Lt=Lt,
                    sel=sel,
                    mat_FM_Lc_com.C=mat_FM_Lc_com.C,
                    mat_FM_Lc_com.Y=mat_FM_Lc_com.Y,
                    mat_FM_Lc_com.V=mat_FM_Lc_com.V,
                    mat_FM_Lc_com.B=mat_FM_Lc_com.B))

    }
  }



  # return results and plot
  class(ret) <- "predict_mod"
  if(plot) plot(ret)
  return(ret)
}

