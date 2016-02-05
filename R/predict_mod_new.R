#' @title Prediction models: Beverton & Holt's yield per recruit model and the Thompson & Bell model
#
#' @description  This is a function predicting catch, yield, biomass and economic values for different
#'    fishing mortality scenarions (in combination with gear changes). It allows to choose if the
#'    calculation should be based upon Beverton and Holt's yield per recruit model or the
#'    more sophisticated Thompson and Bell model.
#'
#' @param param a list consisting of following parameters (not all are required):
#' \itemize{
#'   \item \strong{Linf} or \strong{Winf}: infinite length for investigated species in cm [cm],
#'   \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'   \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \strong{M}: natural mortality,
#'   \item \strong{Z}: total mortality,
#'   \item \strong{FM}: fishing mortality,
#'   \item \strong{a}: ,
#'   \item \strong{b}: ,
#'   \item \strong{Lr} or \strong{tr}: length or age of recuritment;}
#' additional list objects for the Thompson and Bell model:
#'  \itemize{
#'   \item \strong{age} or \strong{midLengths}: midpoints of the length class as vector (length-frequency
#'   data) or ages as vector (age composition data),
#'   \item \strong{meanWeight}: vector with mean weight per length group or age class,
#'   \item \strong{meanValue}: vector with mean value per length group or age class,
#' }
#' @param FM_change vector with ascending fishing mortalities
#' @param Lc_tc_change vector with ascending lengths or ages at first capture (Lc/tc)
#' @param type indicating which model should be applied: \code{"ypr"} for Beverton and Holt's
#'   yield per recruit model and \code{"ThompBell"} for the Thompson and Bell model
#' @param s_list list with selectivity parameters
#' @param stock_size_1 stock size of smallest size class, if NA values are calculated relative to a stock size of 1000 individuals
#' @param unit.time in which time unit the data is provided? "month" or "year"
#' @param plus.group if a value is provided, a plus group is created comprising this size class and all above
#' @param curr.Lc_tc current Lc (length at first capture)
#' @param curr.E current exploitation rate
#'
#' @keywords function prediction
#'
#' @examples
#' \donttest{
#' #______________________________________
#' #       Yiel Per Recruit (YPR)
#' #______________________________________
#' # age structured data
#' # Nemipterus marginatus
#' threadfin <- list(Winf = 286,K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4)
#'
#' # run model
#' predict_mod(threadfin, FM_change = seq(0,6,0.1),
#'    Lc_tc_change = seq(0.2,1,0.2), type = 'ypr')  #where it is maximal  = MSY
#'
#' # Leiognathus spendens (Pauly 1980)
#' ponyfish <- list(Winf = 64, K = 1, t0 = -0.2, M = 1.8, tr = 0.2)
#'
#' # run model
#' predict_mod(ponyfish, Lc_tc_change = c(0.2,0.3,1.0), type = 'ypr')
#'
#' #______________________________________
#' # length structured data
#' # Xiphias gladius (Berkeley and Houde 1980)
#' swordfish <- list(Linf = 309, K = 0.0949, M = 0.18,
#'                   a=0.0003, b=3, Lr = 90, growthFun = "growth_VB")  ## T_Lr , a, b ??? assumed
#'
#' select.list <- list(selecType = 'trawl_ogive', L50 = 120, L75 = 132) ###Lc = 100)
#' #swordfish$midLengths <- seq(60,300,5)
#'
#' # run model
#' predict_mod(param = swordfish, Lc_tc_change = c(100,118,150,180),
#'             s_list = select.list, type = 'ypr', Lmin = 90, Lincr = 8)
#'
#' ####test: E <- seq(0,0.9,0.1) FM <- E * M / (1 - E)
#'
#'
#' #______________________________________
#' #      Thompson and Bell model
#' #______________________________________
#' # with age structured data
#' # load data
#' data(shrimps)
#'
#' # run model
#' predict_mod(shrimps, FM_change = seq(0.1,3,0.1), type = 'ThompBell')
#'
#' # create list with selectivity information
#' select.list <- list(selecType = 'knife_edge',  #or 'gillnet' or 'trawl_ogive'
#'    Lc = 34, tc = 5, selecDist = 'lognormal',    #or 'normal_fixed'
#'    mesh_size = 8.1, mesh_size1 = 9.1, select_p1 = 21.1, select_p2 = 23.8)
#'
#' # add additional parameters to data list
#' shrimps <- c(shrimps,  list(Linf = 50, K = 0.3, t0 = 0.01))
#'
#' # run model
#' predict_mod(shrimps,FM_change = seq(0,3,0.2), Lc_tc_change = seq(24,44,2),
#'    type = 'ThompBell', s_list = select.list, unit.time = 'month')
#'
#' #______________________________________
#' # with length structured data
#' # load data
#' data(hake)
#'
#' # run model
#' predict_mod(hake,FM_change = seq(0.1,3,0.1), type = 'ThompBell')
#' }
#'
#' @details better to treat last group always as a plus group.....
#'    The Thompson and Bell model incorporates an iteration step simulating the stock by means
#'    of the \code{\link{stock_sim}} function. In case changes in gear characteristics -
#'    here measured in terms of Lc or tc, the length or age at first capture, respectively -
#'    should be explored, a list with selectivity information about the gear has to be provided and
#'    the prediction models make use of the selectivity \code{\link{select_ogive}} function.
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
#'      \item FM_change
#'      \item Y_R
#'      \item Y_R.rel
#'      \item B_R
#'      \item B_R.percent
#'      \item Ty
#'      \item Wy
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
#'      \item \strong{tot.C}: total catches for different x factors,
#'      \item \strong{tot.Y}: total yield values for different x factors,
#'      \item \strong{tot.V}: total values for different x factors,
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
#'
#' }
#'
#' @references
#' example 1 : Kuwait (Garcia and van Zalinge 1982)
#'
#' Beverton and Holt 1966
#'
#' Boerema, L.K., and J.A. Gulland, 1973. Stock assessment of the Peruvian anchovy
#' (Engraulis ringens) and management of the fishery. Journal of the Fisheries Board of
#' Canada, 30(12):2226-2235
#'
#' Millar, R.B., and R. Holst, 1997. Estimation of gillnet and hook selectivity using
#' log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3):471-477
#'
#' Sparre, P., and Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' #@export

predict_mod <- function(param, FM_change = NA, Lc_tc_change = NULL, type,  s_list = NA,
                        stock_size_1 = NA, unit.time = 'year', curr.E = NA, curr.Lc_tc = NA,
                        plus.group = NA, Lmin = NA, Lincr = NA){
  res <- param

  # Beverton and Holt's ypr
  #------------
  if(type == "ypr"){
    M <- res$M
    K <- res$K
    t0 <- ifelse(!is.null(res$t0),res$t0,0)
    a <- ifelse(!is.null(res$a),res$a,NA)
    b <- ifelse(!is.null(res$b),res$b,NA)

    if(length(FM_change) == 1 & is.na(FM_change[1])){
      FM_change <- seq(0,10,0.1)
      print(noquote("No fishing mortality (FM_change) was provided, a default range of 0 to 10 is used."))
    }

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #                        Age data                          #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if("Winf" %in% names(res)){
      Winf <- res$Winf
      tr <- res$tr
      tc <- Lc_tc_change

      # yield functions
      # ___________________________________________________
      #Biomass per Recruit
      bpr <- function(x,z){
        S <- exp(-K * (z - t0))
        y <- exp(-M*(z-tr)) * Winf *
          ((1/(x + M)) - ((3*S)/((M + x) + K)) +
             ((3*(S^2))/((M + x)+(2*K))) - ((S^3)/((M + x) + (3*K))))
        return(y)
      }
      # yield per recruit
      ypr <- function(x,z){
        y <- x * (bpr(x,z))
        return(y)
      }
      # realtive yield per recruit
      ypr.rel <- function(x,z){
        y <- ypr(x,z) * (exp(-M *(tr - t0))) / Winf
        return(y)
      }
      # derivative
      derivative <- function(x,z){
        S <- exp(-K * (z - t0))
        C <- ((K*(1-x))/M)
        B <- (S^(M/K)) * (1 - ((3*S)/(1+C)) + ((3*S^2)/(1+2*C)) - ((S^3)/(1+3*C)))
        D <- (-((3*K*S^3)/(M*((3*K*(1-x)/M)+1)^2)) + ((6*K*S^2)/(M*((2*K*(1-x)/M)+1)^2)) -
                ((3*K*S)/(M*((K*(1-x)/M)+1)^2)))
        des <- x * (S^(M/K)) * D  + B
        return(des)
      }
      # ___________________________________________________

      list_tc_runs <- vector("list",length(tc))
      list_Es <- vector("list",length(tc))
      for(i in 1:length(tc)){
        tci <- tc[i]

        Z <- (M + FM_change)
        E <- FM_change/Z

        #Biomass per Recruit
        B_R <- bpr(FM_change, tci)

        #virgin biomass
        Bv_R <- B_R[which(FM_change == 0)]
        #biomass of exploited part of the cohort (biomass of fish older than tc)

        #biomass in percetage of virgin biomass
        B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

        #Yield per Recruit
        #Y_R <- B_R * FM_change
        Y_R <- ypr(FM_change, tci)

        #relative yield per recruit - mostly done with length frequency data (exclusively?)
        #Y_R.rel <- Y_R * (exp(-M *(tr - t0))) / Winf
        Y_R.rel <- ypr.rel(FM_change, tci)

        #mean age in annual yield
        Ty <- (1 / (M+FM_change)) + tci

        #mean length in the annual yield
        #Ly <- Linf * (1 - (((M+FM_change)*S)/((M+FM_change)+K)))

        #mean weight in annual yield
        S <- exp(-K * (tci - t0))
        Wy <- (M+FM_change) * Winf *
          ((1/(FM_change+M)) - ((3*S)/((M+FM_change)+K)) +
             ((3*(S^2))/((M+FM_change)+(2*K))) - ((S^3)/((M+FM_change) + (3*K))))


        results.PBH <- data.frame(FM = FM_change,
                                  Ty = Ty,
                                  Wy = Wy,
                                  Y_R = Y_R,
                                  Y_R.rel = Y_R.rel,
                                  B_R = B_R,
                                  B_R.percent = B_R.percent)


        list_tc_runs[[i]] <- results.PBH

        # First derivative of relative yield per recruit model
        deri <- derivative(E, tci)

        # reference points
        N01 <- which.min(abs(deri - (deri[1] * 0.1)))
        N05 <- which.min(abs(B_R.percent - 50))    ##which.min(abs(deri - (deri[1] * 0.5)))
        Nmax <- which.min(abs(deri))

        df_loop_Es <- data.frame(tc = tci,
                                 F01 = FM_change[N01],
                                 F05 = FM_change[N05],
                                 Fmax = FM_change[Nmax],
                                 E01 = E[N01],
                                 E05 = E[N05],
                                 Emax = E[Nmax])
        list_Es[[i]] <- df_loop_Es
      }

      df_Es <- do.call(rbind,list_Es)

      # current exploitation rate
      curr.F = (M * curr.E)/(1-curr.E)
      df_currents <- data.frame(curr.E = curr.E,
                                curr.F = curr.F,
                                curr.YPR = ypr((curr.F), curr.Lc_tc),
                                curr.YPR.rel = ypr.rel(curr.F, curr.Lc_tc))

      names(list_tc_runs) <- tc
      ret <- c(res,list(
        FM = FM_change,
        tc = tc,
        list_tc_runs = list_tc_runs,
        df_Es = df_Es,
        currents = df_currents))

      class(ret) <- "predict_mod"
      # plot results
      plot(ret)

      return(ret)
    }

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #                       Length data                        #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if("Linf" %in% names(res)){
      Linf <- res$Linf
      Lr <- res$Lr
      Lc <- Lc_tc_change

      # yield functions
      # ___________________________________________________
      # yield function (FM_change, Lci)
      ypr <- function(FM_change,Lci){
        S <- 1 - (Lci/Linf)  # == U  ##(U <- 1 - (Lci/Linf))
        A <- ((Linf - Lci)/(Linf - Lr)) ^ (M/K)
        G <- FM_change + M   # G == Z
        y <- FM_change * A * Winf * ((1/G) - (3*S)/(G + K) + (3*S^2)/(G + 2*K) - (S^3)/(G + 3*K))
        return(y)
      }
      # relative yield function, according to Sparre & Venema and Gayanilo 2006:
      ypr.rel <- function(E, Lci, Lti = NA){
        if(is.na(Lti)) S <- 1 - (Lci/Linf)
        if(!is.na(Lti)) S <- 1 - (Lti/Linf)
        m <- (1-E)/(M/K)    ## == K/Z
        y <- E * (S)^(M/K) * (1 - ((3*S)/(1+m)) + ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))
        # according to Gayanilo 1997 Fisat description (wrong???):
        # y <- E * (S)^(m) * (1 - ((3*S)/(1+m)) + ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))
        return(y)
      }
      # derivative of yield function
      derivative <- function(E, Lci){
        S <- 1 - (Lci/Linf)
        C <- ((K*(1-E))/M)
        B <- (S^(M/K)) * (1 - ((3*S)/(1+C)) + ((3*S^2)/(1+2*C)) - ((S^3)/(1+3*C)))
        D <- (-((3*K*S^3)/(M*((3*K*(1-E)/M)+1)^2)) + ((6*K*S^2)/(M*((2*K*(1-E)/M)+1)^2)) -
                ((3*K*S)/(M*((K*(1-E)/M)+1)^2)))
        des <- E * (S^(M/K)) * D + B
        return(des)
      }
      # relative yield function with selection ogive: i = length classes from Lmin to Lmax
      ypr.rel.sel <- function(expl, pcap, lengths){
        E <- expl
        P <- pcap
        Lt <- lengths
        Y_R.rel.tot.all.classes <- rep(NA,length(E))
        for(Elevels in 1:length(E)){
          Elevel <- E[Elevels]
          # population levels
          # Calculations per size class

          U = 1 - (Lt/Linf)      ###### BIG ASSUMPTION THAT LC = Lt

          # reduction factor per size group
          r <- rep(NA, length(Lt))
          for(x1 in 2:length(Lt)){
            r[x1] <- (U[x1] ^ ((M/K) * (Elevel/(1-Elevel))*P[x1]))  /  (U[x1-1] ^ ((M/K) * (Elevel/(1-Elevel))*P[x1]))
          }

          # G per size group
          G <- rep(NA,length(Lt))
          for(x2 in 1:length(Lt)){
            G[x2] <- prod(r[1:x2], na.rm = TRUE)
          }
          G[1] <- r[1]  # because: rLmin-1 = 1
          G[length(Lt)] <- 0  # because: rLinf = 0

          # Yield per size group
#           Y_R.rel = Elevel * (1 - (Lt/Linf))^(M/K) * (1 - ((3*(1 - (Lt/Linf)))/(1+((1-Elevel)/(M/K)))) +
#                                                          ((3*(1 - (Lt/Linf))^2)/(1+2*((1-Elevel)/(M/K)))) -
#                                                          (((1 - (Lt/Linf))^3)/(1+3*((1-Elevel)/(M/K)))))

          Y_R.rel = Elevel * (1 - (Lt/Linf))^((1-Elevel)/(M/K)) * (1 - ((3*(1 - (Lt/Linf)))/(1+((1-Elevel)/((1-Elevel)/(M/K))))) +
                                                                     ((3*(1 - (Lt/Linf))^2)/(1+2*((1-Elevel)/((1-Elevel)/(M/K))))) -
                                                                     (((1 - (Lt/Linf))^3)/(1+3*((1-Elevel)/((1-Elevel)/(M/K))))))


          Y_R.rel.tot <- rep(NA,length(Lt))
          for(x3 in 2:(length(Lt)-1)){
            Y_R.rel.tot[x3] <- (P[x3]*((Y_R.rel[x3]*G[x3-1]) - (Y_R.rel[x3+1]*G[x3])))
          }

          Y_R.rel.tot.all.classes[Elevels] <- sum(Y_R.rel.tot, na.rm=TRUE)
        }
        return(Y_R.rel.tot.all.classes)
      }
      # relative biomass function with selction ogive
      bpr.rel.sel <- function(expl, pcap, lengths){
        E <- expl
        P <- pcap
        Lt <- lengths
        B_R.rel.tot.all.classes <- rep(NA,length(E))
        for(Elevels in 1:length(E)){
          Elevel <- E[Elevels]
          # population levels
          # Calculations per size class

          U = 1 - (Lt/Linf)      ###### BIG ASSUMPTION THAT LC = Lt

          # reduction factor per size group
          r <- rep(NA, length(Lt))
          for(x1 in 2:length(Lt)){
            r[x1] <- (U[x1] ^ ((M/K) * (Elevel/(1-Elevel))*P[x1]))  /  (U[x1-1] ^ ((M/K) * (Elevel/(1-Elevel))*P[x1]))
          }

          # G per size group
          G <- rep(NA,length(Lt))
          for(x2 in 1:length(Lt)){
            G[x2] <- prod(r[1:x2], na.rm = TRUE)
          }
          G[1] <- r[1]  # because: rLmin-1 = 1
          G[length(Lt)] <- 0  # because: rLinf = 0

          m = (1-Elevel)/(M/K)  #  == K/Z
          mx = m / (1 - Elevel) # == 1/(M/K)
          A <- (1 - ((3*U)/(1+m)) + ((3*U^2)/(1+2*m)) - ((U^3)/(1+3*m)))
          B <- (1 - ((3*U)/(1+mx)) + ((3*U^2)/(1+2*mx)) - ((U^3)/(1+3*mx)))
          B_R.rel <- (1-Elevel) * (A/B)

          B_R.rel.tot <- rep(NA,length(Lt))
          for(x3 in 2:(length(Lt)-1)){
            B_R.rel.tot[x3] <- (P[x3]*((B_R.rel[x3]*G[x3-1]) - (B_R.rel[x3+1]*G[x3])))
          }

          B_R.rel.tot.all.classes[Elevels] <- sum(B_R.rel.tot, na.rm=TRUE)
        }
        return(B_R.rel.tot.all.classes)
      }
      # derivative of selectivity function
      derivative.sel <- function(expl, pcap, lengths){
        E <- expl
        P <- pcap
        Lt <- lengths
        dev.tot.all.classes <- rep(NA,length(E))
        for(Elevels in 1:length(E)){
          Elevel <- E[Elevels]
          # Calculations per size class

          U = 1 - (Lt/Linf)      ###### BIG ASSUMPTION THAT LC = Lt

          # reduction factor per size group
          r <- rep(NA, length(Lt))
          for(x1 in 2:length(Lt)){
            r[x1] <- (U[x1] ^ ((M/K) * (Elevel/(1-Elevel))*P[x1]))  /
              (U[x1-1] ^ ((M/K) * (Elevel/(1-Elevel))*P[x1]))
          }

          # derivative of r
          r.dev <- rep(NA,length(Lt))
          for(x1a in 2:length(Lt)){
            num <- U[x1a] ^ ((M/K) * (Elevel/(1-Elevel))*P[x1a])
            denum <- U[x1a-1] ^ ((M/K) * (Elevel/(1-Elevel))*P[x1a])
            dev.num <-
              (M * P[x1a] *log(U[x1a])) / ((K * U[x1a] ^ ((M * P[x1a] * Elevel) /(K*Elevel - K))) *(Elevel - 1)^2)
            dev.denum <-
              (M * P[x1a] *log(U[x1a-1])) / ((K * U[x1a-1] ^ ((M * P[x1a] * Elevel) /(K*Elevel - K))) *(Elevel - 1)^2)

            r.dev[x1a] <- (dev.num * denum - num * dev.denum) / denum ^2
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
            re.G <- rep(NA,length(x2a))
            for(x2b in 1:x2a){
              if(length(r.dev[x2b] * r[1:x2a][-x2b]) == 0){
                re.G[x2b] <- NA
              }else re.G[x2b] <- r.dev[x2b] * prod(r[1:x2a][-x2b],na.rm=TRUE)
            }
            G.dev[x2a] <- sum(re.G, na.rm=TRUE)
          }

          # Yield per size group
          Y_R.relA = Elevel * (1 - (Lt/Linf))^((1-Elevel)/(M/K)) * (1 - ((3*(1 - (Lt/Linf)))/(1+((1-Elevel)/((1-Elevel)/(M/K))))) +
                                                        ((3*(1 - (Lt/Linf))^2)/(1+2*((1-Elevel)/((1-Elevel)/(M/K))))) -
                                                        (((1 - (Lt/Linf))^3)/(1+3*((1-Elevel)/((1-Elevel)/(M/K))))))
          Y_R.rel.tot <- rep(NA,length(Lt))
          for(x3 in 2:(length(Lt)-1)){
            Y_R.rel.tot[x3] <- (P[x3]*((Y_R.relA[x3]*G[x3-1]) - (Y_R.relA[x3+1]*G[x3])))
          }
          Y_R.rel = Y_R.rel.tot

          # derivative of Y_R.rel per size group
          dev.Y_R.rel <- rep(NA,length(Lt))
          for(x4 in 1:length(Lt)){
            S = 1 - (Lt[x4]/Linf)
            A = (K*Elevel - K)/M
            B = M/K
            C = ((-(3*S)/(B+1)) - ((S^3)/(3*B+1)) + ((3*S^2)/(2*B+1)) + 1)
            dev.Y_R.rel[x4] <- (C / (S^A)) - ((K*log(S)*Elevel*C) / (M*S^A))

            #dev.Y_R.rel[x4] <- derivative(Elevel, Lt[x4])
          }


          # total derivative
          dev.tot <- rep(NA,length(Lt))
          for(x5 in 2:(length(Lt)-1)){
            firstA <- Y_R.rel[x5]*G.dev[x5-1] + dev.Y_R.rel[x5]*G[x5-1]
            secondA <- Y_R.rel[x5+1]*G.dev[x5] + dev.Y_R.rel[x5+1]*G[x5]
            dev.tot[x5] <- P[x5] * (firstA - secondA)
            #dev.tot[x3] <- (P[x3]*((dev.Y_R.rel[x3]*G.dev[x3-1]) - (dev.Y_R.rel[x3+1]*G.dev[x3])))
          }
          dev.tot.all.classes[Elevels] <- sum(dev.tot, na.rm=TRUE)
        }
        return(dev.tot.all.classes)
      }
      # ___________________________________________________

      list_Lc_runs <- vector("list", length(Lc))
      list_Es <- vector("list", length(Lc))
      for(i in 1:length(Lc)){
        Lci <- Lc[i]

        Z <- (M + FM_change)
        E <- FM_change/Z
        # transform Linf in Winf
        Winf <- a * (Linf ^ b)

        # Other selectivity than assumed knife edge
        if(length(s_list) > 1){
          if("midLengths" %in% names(res)){
            Lincr <- res$midLengths[2] - res$midLengths[1]
            Lmin <- res$midLengths[1] - (Lincr/2)
          }
          if(!"midLengths" %in% names(res)){
            Lmin <- Lmin
            Lincr <- Lincr
          } # make Lt fixed except depending on Linf and as detailed as possible, e.g. Lincr = 0.5 (takes a lot time)
          Lt <- seq(Lmin,Linf,Lincr)
          P <- select_ogive(s_list, Lt =  Lt, Lc = Lci)
          Y_R.rel <- ypr.rel.sel(E, P, Lt)

          # convert Y_R.rel to Y_R
          tr = ((log(1-(Lr/Linf)))/K) + t0  # Lr = Linf * (1 - exp(-K * (tr - t0)))
          Y_R = Y_R.rel / (exp(M*(tr - t0)/Winf))  # Y_R.rel  = Y_R * exp(M*(tr - t0)/Winf)
          B_R.rel = bpr.rel.sel(E, P, Lt)  ## NEW

          if(i == length(Lc)) plot(Lt, P, type = 'l', ylab = 'Prob of capture',
                                   main = 'Selectivity function')
        }

        # KNIFE EDGE SELECTION
        #relative yield per recruit - mostly done with length frequency data (exclusively?)
        if(length(s_list) == 1){
          Y_R.rel <- ypr.rel(E, Lci)
          # yield per recruit for length data   # also possbile inputing option: F/K
          Y_R <- ypr(FM_change,Lci)
#           #biomass per recruit for length data?
#           B_R <- Y_R / FM_change
          # relative biomass per recruit for length data?
          B_R.rel <- Y_R.rel / FM_change
          #virgin biomass
          B_R <- Y_R/ FM_change
          Bv_R <- B_R[(which(FM_change == 0)+1)]  ### CHECK: if == 0 than 0 NaN because yield and FM_change is 0
          #biomass of exploited part of the cohort (biomass of fish older than tc)

          #biomass in percetage of virgin biomass
          B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)
        }

        #mean length in the annual yield
        S <- 1 - (Lci/Linf)
        Ly <- Linf * (1 - ((Z*S)/(Z+K)))

        #mean weight in annual yield
        # Wy <- (M+FM_change) * Winf *
        #    ((1/(FM_change+M)) - ((3*S)/((M+FM_change)+K)) +
        #       ((3*(S^2))/((M+FM_change)+(2*K))) - ((S^3)/((M+FM_change) + (3*K))))

        results.PBH <- data.frame(FM = FM_change,
                                  Ly = Ly,
                                  E = E,
                                  Y_R = Y_R,
                                  Y_R.rel = Y_R.rel,
                                  B_R = B_R,
                                  B_R.percent = B_R.percent)
        list_Lc_runs[[i]] <- results.PBH


        # First derivative of relative yield per recruit model
        if(length(s_list) == 1) deri <- derivative(E, Lci)
        if(length(s_list) > 1) deri <- derivative.sel(E, P, Lt)

        # reference points
        N01 <- which.min(abs(deri - (deri[1] * 0.1)))
        N05 <- which.min(abs(deri - (deri[1] * 0.5)))   ##which.min(abs(B_R.percent - 50)) # wrong!!!
        Nmax <- which.min(abs(deri))


        df_loop_Es <- data.frame(Lc = Lci,
                                 F01 = FM_change[N01],
                                 F05 = FM_change[N05],
                                 Fmax = FM_change[Nmax],
                                 E01 = E[N01],
                                 E05 = E[N05],
                                 Emax = E[Nmax]
        )
        list_Es[[i]] <- df_loop_Es
      }

      df_Es <- do.call(rbind,list_Es)

      # current exploitation rate
      curr.F = (M * curr.E)/(1-curr.E)
      df_currents <- data.frame(curr.E = curr.E,
                                curr.F = curr.F,
                                curr.YPR = ypr((curr.F), curr.Lc_tc),
                                curr.YPR.rel = ypr.rel(curr.E, curr.Lc_tc))

      names(list_Lc_runs) <- Lc
      ret <- c(res,list(FM = FM_change,
                        Lc = Lc,
                        list_Lc_runs = list_Lc_runs,
                        df_Es = df_Es,
                        currents = df_currents))
      class(ret) <- "predict_mod"

      # plot results
      #####plot(ret)

      return(ret)
    }
  }

  #------------


  # Thompson and Bell model
  #------------
  if(type == "ThompBell"){
    meanWeight <- res$meanWeight
    meanValue <- res$meanValue

    #mortalities
    FM <- res$FM
    if(!is.null(res$M)){
      nM <- res$M
      Z <- FM + nM
    }else{
      Z <- res$Z
      nM <- mean(Z - FM,na.rm=T)
    }

    Linf <- res$Linf
    K <- res$K
    t0 <- ifelse("t0" %in% names(res),res$t0,0)

    #   FM <- res$FM
    #   Z <- res$Z
    #   #natural Mortality
    #   nM <- mean(Z - FM,na.rm=T)

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #                        Age data                          #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if('age' %in% names(res)) classes <- as.character(res$age)
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #                       Length data                        #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if('midLengths' %in% names(res)) classes <- as.character(res$midLengths)

    #HHHHHHHHHHHH#
    #     ALL    #
    #HHHHHHHHHHHH#

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    # Only FM change provided without Lc_tc change
    if(is.null(Lc_tc_change)){
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
      class(ret) <- "predict_mod"

      # plot results
      plot(ret)

      return(ret)
    }


    # FM and Lc_tc change provided
    if(!is.null(Lc_tc_change)){
      # instead of s_list the outcome of one of the other select functions?
      #as a option or put values per hand

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

      pred.FM_Lc_com_res_loopC_list <- list()
      pred.FM_Lc_com_res_loopY_list <- list()
      pred.FM_Lc_com_res_loopB_list <- list()
      pred.FM_Lc_com_res_loopV_list <- list()

      for(x21 in 1:length(FM_Lc_com_mat.list)){  #loop for length of list == Lc changes
        mati <- FM_Lc_com_mat.list[[x21]]

        pred.FM_Lc_com_res_loop1_list <- list()
        for(x22 in 1:dim(mati)[2]){

          param.loop$FM <- mati[,x22]
          param.loop$Z <- mati[,x22] + nM
          res2 <- stock_sim(param.loop, unit.time,
                            stock_size_1,plus.group=plus.group)
          pred.FM_Lc_com_res_loop1_list[[x22]] <- res2$totals
        }
        prev_mat <- do.call(rbind, pred.FM_Lc_com_res_loop1_list)
        prev_matC <- prev_mat[,'tot.C']
        prev_matY <- prev_mat[,'tot.Y']
        prev_matB <- prev_mat[,'meanB']
        prev_matV <- prev_mat[,'tot.V']

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
                    Lt=Lt,sel=sel,
                    mat_FM_Lc_com.C=mat_FM_Lc_com.C,
                    mat_FM_Lc_com.Y=mat_FM_Lc_com.Y,
                    mat_FM_Lc_com.V=mat_FM_Lc_com.V,
                    mat_FM_Lc_com.B=mat_FM_Lc_com.B))
      class(ret) <- "predict_mod"

      # plot results
      plot(ret)

      return(ret)
    }
  }
  #------------

}





## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?

