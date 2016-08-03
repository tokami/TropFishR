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
#' \dontrun{
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
#' }
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

predict_mod <- function(param, type, FM_change = NA, E_change = NA, Lc_change = NULL,
                        tc_change = NULL,
                        s_list = NA,
                        stock_size_1 = NA, age_unit = 'year', curr.E = NA,
                        curr.Lc_tc = NA,
                        plus_group = NA, Lmin = NA, Lincr = NA, plot = FALSE){
  res <- param

  # Beverton and Holt's ypr
  if(type == "ypr"){
    M <- res$M
    K <- res$K
    t0 <- ifelse(!is.null(res$t0),res$t0,0)
    a <- res$a  # might be NULL
    b <- res$b  # might be NULL
    Winf <- res$Winf  # might be NULL
    Linf <- ifelse(!is.null(res$Linf),res$Linf, exp(log(Winf/a)/b))
    # REALLY ? maybe without Linf: then message that Winf has to exist
    if(is.null(Linf) | is.na(Linf)) stop("Either Linf or Winf with a and b has to be provided!")
    if(is.null(Winf)) Winf <-  a * Linf ^ b###exp((log(Linf) - a)/b) # might still be NULL
    # or              Winf <- exp(log(Linf-a)/b)



    if(length(FM_change) == 1 & is.na(FM_change[1]) & length(E_change) == 1 & is.na(E_change[1])){
      FM_change <- seq(0,10,0.1)
      print(noquote("No fishing mortality (FM_change) or exploitation rate (E_change) was provided, a default range for fishing mortality of 0 to 10 is used."))
    }

    # transfer E_change into F_change if provided
    if(length(FM_change) == 1 & is.na(FM_change[1]) & length(E_change) != 1 & !is.na(E_change[1])){
      FM_change <- (E_change * M) / (1 - E_change)
      FM_change[FM_change == Inf] <- (0.9999 * M) / (1 - 0.9999)
    }

    # Recruitment  - knife edge
    tr <- res$tr   # might be NULL
    Lr <- res$Lr   # might be NULL
    if(is.null(tr) & is.null(Lr)) stop("Either the age or the length at recruitment (tr or Lr) has to be provided in param!")
    if(is.null(tr)) tr <- VBGF(L=Lr,param = list(Linf=Linf,K=K,t0=t0)) # VBGF(L=Lr,Linf=Linf,K=K,t0=t0)
    if(is.null(Lr)) Lr <- VBGF(t=tr,param = list(Linf=Linf,K=K,t0=t0)) # VBGF(t=tr,Linf=Linf,K=K,t0=t0)


    # Selectivity - knife edge or with selctivtiy ogive
    tc <- res$tc   # might be NULL
    Lc <- res$Lc   # might be NULL
    if(is.null(tc) & is.null(Lc)){
      if("L50" %in% s_list) Lc <- s_list$L50
      if("Lc" %in% s_list) Lc <- s_list$Lc
      if(!("Lc" %in% s_list) & !("L50" %in% s_list))stop("Either the age or the length at first capture (tc or Lc) has to be provided in param! \n Or provide a Lc value in s_list!")
    }
    if(is.null(tc)) tc <- VBGF(L=Lc, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(L=Lc,Linf=Linf,K=K,t0=t0)
    if(is.null(Lc)) Lc <- VBGF(t=tc, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(t=tc,Linf=Linf,K=K,t0=t0)
    if(is.null(tc_change) & !is.null(Lc_change)) tc_change <- VBGF(L=Lc_change, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(L=Lc_change,Linf=Linf,K=K,t0=t0)
    if(is.null(Lc_change) & !is.null(tc_change)) Lc_change <- VBGF(t=tc_change, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(t=tc_change,Linf=Linf,K=K,t0=t0)
    tc <- c(tc,tc_change)
    Lc <- c(Lc,Lc_change)

    if(length(s_list) > 1){
      selecType <- s_list$selecType
    }else{
      selecType <- "knife_edge"
    }


    # HEART
    list_Lc_runs <- vector("list", length(Lc))
    list_Es <- vector("list", length(Lc))

    # show progress bar only if the loop has more than 1 runs
    nlk <- length(Lc)
    if(nlk > 1){
      pb <- txtProgressBar(min=1, max=nlk, style=3)
      counter <- 1
    }

    for(i in 1:length(Lc)){

      Lci <- Lc[i]
      tci <- tc[i]

      Z <- (M + FM_change)
      E <- FM_change/Z

      # KNIFE EDGE
      if(length(s_list) == 1 | selecType == "knife_edge"){
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
          if(is.na(Lmin) | is.na(Lincr)) writeLines("No midpoints of length classes are provided. This can be done using Lmin and Lincr or by a \n midLengths element in param. A standard range from 1cm to Linf * 0.98 by 2cm is assumed")
          Lmin <- ifelse(is.na(Lmin), 1, Lmin)
          Lincr <- ifelse(is.na(Lincr), 2, Lincr)
          # mid length vector
          Lt <- seq(Lmin,(Linf*0.98),Lincr)
        } # make Lt fixed except depending on Linf and as detailed as possible, e.g. Lincr = 0.5 (takes a long time)

        # selectivity
        P <- select_ogive(s_list, Lt =  Lt, Lc = Lci)

        input <- list(Linf=Linf,
                      Winf = Winf,
                      K = K,
                      M = M,
                      t0 = t0,
                      tr = tr,
                      tc = tci)
        output <- ypr_sel(input, FM_change, Lt, P)

        # relative yield and biomass per recruit
        B_R.rel <- output$rbr
        Y_R.rel <- output$ryr
        # derivative
        deri <- output$derivative

        #test
        Y_R <- Y_R.rel * Winf * exp(M * (tr - t0))
        B_R <- B_R.rel * Winf * exp(M * (tr - t0))

        # biased because only prints P for largest Lc value
        #if(i == length(Lc)) plot(Lt, P, type = 'l', ylab = 'Prob of capture',
        #                         main = 'Selectivity function')
      }


      # virgin biomass
      if(0 %in% FM_change){
        Bv_R <- B_R[FM_change == 0]
      }else{
        Bv_R <- B_R[FM_change == min(FM_change,na.rm = TRUE)]
        writeLines(paste0("Biomass was not estimated for a fishing mortality (FM) of 0, thus the virgin biomass corresponds to a FM of ",min(FM_change,na.rm = TRUE)))
      }

      #biomass in percetage of virgin biomass
      B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

      #mean age in annual yield
      Ty <- (1 / Z) + tci

      #mean length in the annual yield
      S <- exp(-K * (tci - t0))         # the same:    S <- 1 - (Lci/Linf)
      Ly <- Linf * (1 - ((Z*S)/(Z+K)))

      #mean weight in annual yield
      Wy <- (Z) * Winf *
        ((1/Z) - ((3*S)/(Z+K)) +
           ((3*(S^2))/(Z+(2*K))) - ((S^3)/(Z + (3*K))))


      results.PBH <- data.frame(FM = FM_change,
                                E = E,
                                Ty = Ty,
                                Ly = Ly,
                                Wy = Wy,
                                Y_R.rel = Y_R.rel,
                                B_R.rel = B_R.rel)
      # WHY NECESSARY???
      if(length(Y_R) > 0) results.PBH$Y_R = Y_R
      if(length(B_R) > 0) results.PBH$B_R = B_R
      if(length(B_R.percent) > 0) results.PBH$B_R.percent = B_R.percent

      list_Lc_runs[[i]] <- results.PBH


      # reference points
      Nmax <- which.min(abs(deri))
      deri_pot <- deri[1:Nmax]
      N01 <- which.min(abs(deri_pot - (deri[1] * 0.1)))
      Nmsy <- which.max(Y_R.rel)
      N05 <- which.min(abs(B_R.percent - 50))  #which.min(abs(deri - (deri[1] * 0.5)))

      df_loop_Es <- data.frame(Lc = Lci,
                               tc = tci,
                               F01 = FM_change[N01],
                               Fmsy = FM_change[Nmsy])
      if(length(B_R.percent) > 0) df_loop_Es$F05 <- FM_change[N05]   # WHY NECESSARY????
      df_loop_Es$Fmax <- FM_change[Nmax]
      df_loop_Es$E01 <- E[N01]
      df_loop_Es$Emsy <- E[Nmsy]
      if(length(B_R.percent) > 0) df_loop_Es$E05 <- E[N05]    # WHY NECESSARY????
      df_loop_Es$Emax <- E[Nmax]

      list_Es[[i]] <- df_loop_Es

      # update counter and progress bar
      if(nlk > 1){
        setTxtProgressBar(pb, counter)
        counter <- counter + 1
      }
    }

    df_Es <- do.call(rbind,list_Es)

    names(list_Lc_runs) <- paste0("Lc_",Lc)   # names(list_tc_runs) <- tc
    ret <- c(res,list(FM = FM_change,
                      Lc = Lc,           #   tc = tc,
                      list_Lc_runs = list_Lc_runs,   #   list_tc_runs = list_tc_runs,
                      df_Es = df_Es))   #   df_Es = df_Es,


    if(!is.na(curr.E) & !is.na(curr.Lc_tc)){
      # current exploitation rate
      curr.F = (M * curr.E)/(1-curr.E)  # curr.F <- (M * curr.E)/(1-curr.E)
      df_currents <- data.frame(curr.E = curr.E,
                                curr.F = curr.F,
                                curr.YPR = ypr(curr.F, curr.Lc_tc),           #, type = "length"),           # curr.YPR = ypr(curr.F, curr.Lc_tc, type = "age"),
                                curr.YPR.rel = ypr.rel(curr.F, curr.Lc_tc),   #, type = "length"),   # curr.YPR.rel = ypr.rel(curr.F, curr.Lc_tc, type = "age"),
                                curr.BPR = bpr(curr.F, curr.Lc_tc),           #, type = "length"),           # curr.BPR = bpr(curr.F, curr.Lc_tc, type = "age"),
                                curr.BPR.rel = bpr.rel(curr.F, curr.Lc_tc))   #, type = "length"))   # curr.BPR.rel = bpr.rel(curr.F, curr.Lc_tc, type = "age"))
      ret$currents <- df_currents
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

