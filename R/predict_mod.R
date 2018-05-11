#' @title Prediction models
#'
#' @description  This function applies Beverton & Holt's yield per recruit model
#'    as well as the Thompson & Bell model. These models predict catch, yield, biomass
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
#' @param LW_unit a character indicating the unit of the LWa parameter, either "kg" for
#'    kg/cm3 or "g" for g/cm3. Default is "g".
#' @param plot logical; if TRUE results are displayed graphically
#' @param mark logical; if value of choosen points should be displayed in graph (default: TRUE)
#' @param hide.progressbar logical; should progressbar be displayed or hidden? (Default: FALSE)
#' @param boot an object of class 'lfqBoot'
#' @param natMort name of column with natural mortalites for bootstrapping application of VPA
#' @param yearSel optional character for bootsrapping use of this function
#'    specifying a year to be subsetted if LFQ data covers multiple years
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
#'    tc_change = seq(0.2,1,0.2), type = 'ypr')  #where it is maximal  = MAX
#'
#' # Leiognathus spendens (Pauly, 1980)
#' ponyfish <- list(Winf = 64, K = 1, t0 = -0.2, M = 1.8, tr = 0.2)
#'
#' predict_mod(ponyfish, tc_change = c(0.2,0.3,1.0), type = 'ypr', plot=TRUE)
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
#' output <- predict_mod(param = swordfish, Lc_change = c(100,118,150,180),
#'             s_list = select.list, type = 'ypr', Lmin = 90, Lincr = 8)
#' plot(output)
#'
#' data(hake)
#' hake$Lr <- 35
#' select.list <- list(selecType = 'trawl_ogive', L50 = 50, L75 = 54)
#' output <- predict_mod(param = hake, FM_change = seq(0,3,0.05),
#'                       Lc_change = seq(30,70,1), s_list = select.list,
#'                       type = 'ypr', plot = FALSE, curr.Lc = 50, curr.E = 0.73)
#' plot(output, type = "Isopleth", xaxis1 = "FM", yaxis1 = "Y_R.rel", mark = TRUE)
#'
#' output <- predict_mod(param = hake, E_change = seq(0,1,0.1),
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
#' output <- predict_mod(param = shrimps, FM_change = seq(0.1,20,0.1),
#'      type = "ThompBell", age_unit = "month", plot = TRUE)
#'
#' #______________________________________
#' # with length structured data
#' data(hake)
#' par(mar = c(5, 4, 4, 7))
#' predict_mod(param = hake,FM_change = seq(0.1,3,0.05),
#'      type = 'ThompBell', plot = TRUE)
#'
#' # create list with selectivity information
#' select.list <- list(selecType = 'trawl_ogive', L50 = 50, L75 = 55)
#'
#' output <- predict_mod(param = hake, FM_change = seq(0,2,0.1),
#'      Lc_change = seq(20,70,1),
#'      curr.E = 0.4, curr.Lc = 50,
#'      type = 'ThompBell', s_list = select.list)
#' plot(output, xaxis1 = "FM", yaxis_iso = "Lc", yaxis1 = "B_R", mark = TRUE)
#'
#' #______________________________________
#' # YPR with bootstrapping
#' # coming soon!
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
#' @importFrom ks kde
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

predict_mod <- function(
  param, type, FM_change = NA,
  E_change = NA,
  FM_relative = FALSE,
  Lc_change = NULL,
  tc_change = NULL,
  s_list = NA,
  stock_size_1 = NA, age_unit = 'year', curr.E = NA,
  curr.Lc = NA,
  plus_group = FALSE, Lmin = NA, Lincr = NA,
  LW_unit = "g",
  plot = FALSE, mark = TRUE,
  hide.progressbar = FALSE,
  boot = NULL,
  natMort = NULL,
  yearSel = NA
){

    ## VPA with bootstrapping ELEFAN results
    if(!is.null(boot) & class(boot) == "lfqBoot"){

        bootOut <- boot
        bootRaw <- boot$bootRaw


        if (any(is.null(bootRaw$Linf), is.null(bootRaw$K)))
            stop("YPR with boot requires a boot object with columns Linf and K")

        if(!("a" %in% names(param)) | !("b" %in% names(param)))
            stop("YPR requires information about the length-weight relationship. Please provide 'a' and 'b' estimates in param.")
        a <- param$a
        b <- param$b
        ## this makes sure that a is always in kg/cm3 = and therefore meanBodyWeight is always in kg
        if(LW_unit == "g"){
            a <- a / 1000
        }


        ## add possibility to provide selecitvity information with s_list
        if(!("L50" %in% names(bootRaw)) &
           !("FMvecVPA" %in% names(boot))){
            stop("YPR requires information about the length at recruitment and length at first capture. Please provide 'Lr' and 'Lc' estimates in param.")
        }


        ## resample data sets
        lfqAll <- lfqResample(param, boot)
        ##  set.seed(boot$seed[bi])


        F01 <- vector("numeric",nrow(bootRaw))
        Fmax <- vector("numeric",nrow(bootRaw))
        F05 <- vector("numeric",nrow(bootRaw))
        for(bi in 1:nrow(bootRaw)){

            lfqTemp <- lfqAll[[bi]]

            ## subset data for specific years (optional) due to seed values not possible before resampling
            if(!is.na(yearSel)){
                yearSel <- as.character(yearSel)
                ## warning if year not in dates
                dateYears <- format(lfqTemp$dates, "%Y")
                if(all(yearSel %in% dateYears == FALSE)) stop(paste0("The selected year ", yearSel, " is not in the LFQ data!"))

                lfqTemp$catch <- lfqTemp$catch[,which(format(lfqTemp$dates,"%Y") %in% yearSel)]
                lfqTemp$dates <- lfqTemp$dates[format(lfqTemp$dates,"%Y") %in% yearSel]
            }

            lfqLoop <- lfqModify(lfqTemp, vectorise_catch = TRUE)


            Linf <- bootRaw$Linf[bi]
            K <- bootRaw$K[bi]
            t_anchor <- bootRaw$t_anchor[bi]
            C <- ifelse("C" %in% names(bootRaw),bootRaw$C[bi],0)
            ts <- ifelse("ts" %in% names(bootRaw),bootRaw$ts[bi],0)
            t0 <- 0

            if(!(natMort %in% names(bootRaw))) stop("Please provide a natural mortality estimate 'M' in the boot object.")
            M <- bootRaw[bi,which(colnames(bootRaw) == natMort)] ## vector with Ms not yet implemented for boot YPR

            if(!("FM" %in% names(bootRaw)))
                stop("Please provide a fishing mortality estimate 'FM' in the boot object.")
            FM <- bootRaw$FM[bi]
            if(is.na(bootRaw$FM[bi]) | is.nan(bootRaw$FM[bi]) |
               is.null(bootRaw$FM[bi]) | bootRaw$FM[bi] == Inf |
               bootRaw$FM[bi] == -Inf){
                    FM <- 0.1  ## hack but might be NaN
            }


            if(!("Z" %in% names(bootRaw)))
                stop("Please provide a total mortality estimate 'Z' in the boot object.")
            Z <- bootRaw$Z[bi]
            if(is.na(bootRaw$Z[bi]) | is.nan(bootRaw$Z[bi]) |
               is.null(bootRaw$Z[bi]) | bootRaw$Z[bi] == Inf |
               bootRaw$Z[bi] == -Inf){
                    Z <- 0.9  ## hack but might be NaN
            }


            catch <- lfqLoop$catch


            if(FM_relative) stop(noquote("ypr does not work with relative changes in FM, please provide absolute values."))


            Winf <- a * (Linf ^ b)


            if(length(FM_change) == 1 & is.na(FM_change[1]) &
               length(E_change) == 1 & is.na(E_change[1])){
              FM_change <- seq(0,10,0.1)
              print(noquote("No fishing mortality (FM_change) or exploitation rate (E_change) was provided, a default range for fishing mortality of 0 to 10 is used."))
            }

            if(length(FM_change) == 1 & is.na(FM_change[1]) &
               length(E_change) != 1 & !is.na(E_change[1])){
              E_change <- E_change[E_change <= 0.9]
              FM_change <- (E_change * M) / (1 - E_change)
            }


            ## Selectivity - knife edge or with selctivtiy ogive
            tc <- param$tc   # might be NULL
            Lc <- param$Lc   # might be NULL
            if("L50" %in% names(bootRaw) & !("Lc" %in% names(param))){
                Lc <- bootRaw$L50[bi]
                if(is.na(bootRaw$L50[bi]) | is.nan(bootRaw$L50[bi]) |
                   is.null(bootRaw$L50[bi]) | bootRaw$L50[bi] == Inf |
                   bootRaw$L50[bi] == -Inf | bootRaw$L50[bi] < 0 |
               bootRaw$L75[bi] >= bootRaw$Linf[bi]){
                    Lc <- 5  ## hack but might be NaN
                }
            }
            if(is.null(tc) & is.null(Lc)){
              if("L50" %in% s_list) Lc <- s_list$L50
              if("Lc" %in% s_list) Lc <- s_list$Lc
              #if(!("Lc" %in% s_list) & !("L50" %in% s_list))stop("Either the age or the length at first capture (tc or Lc) has to be provided in param! \n Or provide a Lc value in s_list!")
            }
            if(!is.null(Linf)){
                if(is.null(tc) & !is.null(Lc)) tc <- VBGF(L=Lc, param = list(Linf=Linf,K=K,t0=t0,ts=ts,C=C))
                ## VBGF(L=Lc,Linf=Linf,K=K,t0=t0)
                if(is.null(Lc) & !is.null(tc)) Lc <- VBGF(t=tc, param = list(Linf=Linf,K=K,t0=t0,ts=ts,C=C))
                ## VBGF(t=tc,Linf=Linf,K=K,t0=t0)
                if(is.null(tc_change) & !is.null(Lc_change))
                    tc_change <- VBGF(L=Lc_change,
                                      param = list(Linf=Linf,K=K,t0=t0,ts=ts,C=C))
                ## VBGF(L=Lc_change,Linf=Linf,K=K,t0=t0)
                if(is.null(Lc_change) & !is.null(tc_change))
                    Lc_change <- VBGF(t=tc_change,
                                      param = list(Linf=Linf,K=K,t0=t0,ts=ts,C=C))
                ## VBGF(t=tc_change,Linf=Linf,K=K,t0=t0)
            }
            tc <- c(tc,tc_change)
            Lc <- c(Lc,Lc_change)


            ## THOMPSON BELL MODEL
            if(type == "ThompBell"){

                ## length based
                classes <- as.character(param$midLengths)

                ## create column without plus group (sign) if present
                classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
                classes.num <- as.numeric(classes.num[,1])
                Lt <- classes.num  # if age in names(res) Lt here is age because classes.num = age

                ## selectivity information  (add possibility to provide externally estimated L50 and L75 in list)
                ## checking for L50/L75 estimation in catch curve
                if("FMvecVPA" %in% names(boot)){
                    FM <- boot$FMvecVPA[[bi]]
                }else if(is.na(s_list) &
                         "L50" %in% names(bootRaw) &
                         "L75" %in% names(bootRaw)){
                    s_list <- list(selecType = "trawl_ogive",
                                   L50 = bootRaw$L50[bi],
                                   L75 = bootRaw$L75[bi])
                    selecType = s_list$selecType
                    if(is.na(bootRaw$L50[bi]) | is.nan(bootRaw$L50[bi]) |
                       is.null(bootRaw$L50[bi]) | bootRaw$L50[bi] == Inf |
                       bootRaw$L50[bi] == -Inf | bootRaw$L50[bi] < 0 |
                       bootRaw$L75[bi] >= bootRaw$Linf[bi]){
                        s_list$L50 <- 5  ## hack
                        s_list$L75 <- 6
                    }
                    sel <- select_ogive(s_list, Lt = Lt)
                    FM <- bootRaw$FM[bi] * sel
                }else{
                    stop("You need to run VPA or catchCurve with calc_ogive to have selectivity information")
                }

                ## add parameters to lfq data
                lfqLoop$Linf <- Linf
                lfqLoop$K <- K
                lfqLoop$t_anchor <- t_anchor
                lfqLoop$C <- C
                lfqLoop$ts <- ts
                lfqLoop$t0 <- t0
                lfqLoop$a <- a
                lfqLoop$b <- b

                #prediction based on f_change
                pred_mat <- as.matrix(FM/max(FM, na.rm = TRUE)) %*% FM_change
                pred_res_list <- vector("list", length(FM_change))
                for(x7 in 1:length(FM_change)){
                  lfqLoop$Z <- pred_mat[,x7] + M
                  lfqLoop$FM <- pred_mat[,x7]
                  resL <- stock_sim(lfqLoop, age_unit = age_unit,
                                    stock_size_1 = stock_size_1, plus_group = plus_group)
                  pred_res_list[[x7]] <- resL$totals
                }

                pred_res_df <- do.call(rbind, pred_res_list)
                pred_res_df$FM_change <- FM_change
                pred_res_df$E_change <- E_change

                res2 <- pred_res_df
                res3 <- c(res,res2)

                
                ## reference points
                ## Fmax
                Nmax <- which.max(pred_res_df$totY)                
                ## F01 (proxy for Fmsy in ICES)
                slopeOrg <- (pred_res_df$totY[2] - pred_res_df$totY[1]) / (FM_change[2] - FM_change[1])
                slope01 <- round(0.1*slopeOrg, 2)
                slopes <- rep(NA, length(FM_change))
                slopes[1] <- slopeOrg
                for(i in 3:length(FM_change)){
                    slopes[i-1] <- round((pred_res_df$totY[i] - pred_res_df$totY[i-1]) /
                        (FM_change[i] - FM_change[i-1]),2)
                }
                dif <- abs(slopes - slope01)
                dif[is.na(dif)] <- 1e+11
                difpot <- dif[1:Nmax]
                N01 <- which.min(difpot)                
                ## F05
                Bper <- rep(NA,length(pred_res_df$meanB))
                Bper[1] <- 100
                for(ix in 2:length(Bper)){
                  Bper[ix] <- pred_res_df$meanB[ix]/pred_res_df$meanB[1] * 100
                }
                N05 <- which.min(abs(Bper - 50))

                ## ref level
                F01[bi] <- FM_change[N01]
                Fmax[bi] <- FM_change[Nmax]
                F05[bi] <- FM_change[N05]


                ## if FMvecVPA not used and selectivity parameters could not be estimated
                if("L50" %in% names(bootRaw)){
                    if(is.na(bootRaw$L50[bi]) | is.nan(bootRaw$L50[bi]) |
                       is.null(bootRaw$L50[bi]) | bootRaw$L50[bi] == Inf |
                   bootRaw$L50[bi] == -Inf | bootRaw$L50[bi] < 0 |
                   bootRaw$L75[bi] >= bootRaw$Linf[bi]){
                    F01[bi] <- NaN  ## hack
                    Fmax[bi] <- NaN
                    F05[bi] <- NaN
                    }
                }
                if(is.na(bootRaw$Z[bi]) | is.nan(bootRaw$Z[bi]) |
                   is.null(bootRaw$Z[bi]) | bootRaw$Z[bi] == Inf |
                   bootRaw$Z[bi] == -Inf){
                    F01[bi] <- NaN  ## hack
                    Fmax[bi] <- NaN
                    F05[bi] <- NaN
                }
            }



            ## BEVERTON HOLT YPR MODEL
            if(type == "ypr"){

                if(!("Lr" %in% names(param)) | (!("Lc" %in% names(param)) & !("L50" %in% names(bootRaw))))
                    stop("YPR requires information about the length at recruitment and length at first capture. Please provide 'Lr' and 'Lc' estimates in param.")
                Lr <- param$Lr
                if("Lc" %in% names(param) & !("L50" %in% names(bootRaw))){
                    Lc <- param$Lc
                }



                # Recruitment  - knife edge
                tr <- param$tr   # might be NULL
                Lr <- param$Lr   # might be NULL
                if(is.null(tr) & is.null(Lr)) stop("Either the age or the length at recruitment (tr or Lr) has to be provided in param!")
                if(!is.null(Linf)){
                    if(is.null(tr)) tr <- VBGF(L=Lr,
                                               param = list(Linf=Linf,K=K,t0=t0,ts=ts,C=C))
                    ## VBGF(L=Lr,Linf=Linf,K=K,t0=t0)
                    if(is.null(Lr)) Lr <- VBGF(t=tr,
                                               param = list(Linf=Linf,K=K,t0=t0,ts=ts,C=C))
                    ## VBGF(t=tr,Linf=Linf,K=K,t0=t0)
                }

                ## checking for L50/L75 estimation in catch curve
                if(is.na(s_list) & "L50" %in% names(bootRaw) & "L75" %in% names(bootRaw)){
                    s_list <- list(selecType = "trawl_ogive",
                                   L50 = bootRaw$L50[bi],
                                   L75 = bootRaw$L75[bi])
                    selecType = s_list$selecType
                    if(is.na(bootRaw$L50[bi]) | is.nan(bootRaw$L50[bi]) |
                       is.null(bootRaw$L50[bi]) | bootRaw$L50[bi] == Inf |
                       bootRaw$L50[bi] == -Inf | bootRaw$L50[bi] < 0 |
                   bootRaw$L75[bi] >= bootRaw$Linf[bi]){
                        s_list$L50 <- 5  ## hack
                        s_list$L75 <- 6
                    }
                }else{
                    if(length(s_list) > 1){
                        selecType <- s_list$selecType
                    }else{
                        selecType <- "knife_edge"
                    }
                }


                # HEART
                list_Lc_runs <- vector("list", length(Lc))
                list_Es <- vector("list", length(Lc))

                if(is.null(Lc)) Lc_tc <- tc else Lc_tc <- Lc
                for(i in 1:length(Lc_tc)){

                  Lci <- Lc[i]
                  tci <- tc[i]

                  Z <- (M + FM_change)
                  E <- FM_change/Z

                  # KNIFE EDGE
                  if(length(s_list) == 1 ){#| selecType == "knife_edge"){
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
                  if(length(s_list) > 1 ){#& selecType != "knife_edge"){
                    if("midLengths" %in% names(param)){
                      classes <- as.character(param$midLengths)
                      # create column without plus group (sign) if present
                      classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
                      classes.num <- as.numeric(classes.num[,1])
                      Lt <- classes.num
                      # Lincr <- param$midLengths[2] - param$midLengths[1]
                      # Lmin <- param$midLengths[1] - (Lincr/2)
                    }
                    if(!"midLengths" %in% names(param)){
                      if(is.na(Lmin) | is.na(Lincr)) writeLines("No midpoints of length classes are provided. This can be done using Lmin and Lincr or by a \n midLengths element in param. A standard range from 1cm to Linf * 0.98 by 2cm is assumed")
                      Lmin <- ifelse(is.na(Lmin), 1, Lmin)
                      Lincr <- ifelse(is.na(Lincr), 2, Lincr)
                      # mid length vector
                      Lt <- seq(Lmin,(Linf*0.98),Lincr)
                    } # make Lt fixed except depending on Linf and as detailed as possible, e.g. Lincr = 0.5 (takes a long time)

                    # selectivity
                    P <- select_ogive(s_list, Lt =  Lt, Lc = Lci)

                    input <- list(Linf = ifelse(length(Linf) > 0, Linf, NA),
                                  Winf = ifelse(length(Winf) > 0, Winf, NA),
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

                  ##biomass in percetage of virgin biomass
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
                                            E = E)
                  if(length(Ty) > 0) results.PBH$Ty <- Ty
                  if(length(Ly) > 0) results.PBH$Ly <- Ly
                  if(length(Wy) > 0) results.PBH$Wy <- Wy

                  results.PBH$Y_R.rel <- Y_R.rel
                  results.PBH$B_R.rel <- B_R.rel

                  # WHY NECESSARY???
                  if(length(Y_R) > 0) results.PBH$Y_R = Y_R
                  if(length(B_R) > 0) results.PBH$B_R = B_R
                  if(length(B_R.percent) > 0) results.PBH$B_R.percent = B_R.percent

                  list_Lc_runs[[i]] <- results.PBH

                  # reference points
                  Nmax <- which.max(Y_R.rel)  #  should be the same as which.min(abs(deri)) which is also labelled Nmax
                  deri_pot <- deri[1:Nmax]
                  N01 <- which.min(abs(deri_pot - (deri[1] * 0.1)))
                    N05 <- which.min(abs(B_R.percent - 50))  #which.min(abs(deri - (deri[1] * 0.5)))


                  df_loop_Es <- data.frame(Lc = ifelse(!is.null(Lci),Lci,NA),
                                           tc = ifelse(!is.null(tci),tci,NA),
                                           F01 = FM_change[N01],
                                           Fmax = FM_change[Nmax])
                  if(length(B_R.percent) > 0) df_loop_Es$F05 <- FM_change[N05]
                  # df_loop_Es$Fmax <- FM_change[Nmax]
                  df_loop_Es$E01 <- E[N01]
                  df_loop_Es$Emax <- E[Nmax]
                  if(length(B_R.percent) > 0) df_loop_Es$E05 <- E[N05]
                  # df_loop_Es$Emax <- E[Nmax]

                  list_Es[[i]] <- df_loop_Es
                }

                curr.E <- FM/Z
                curr.Lc <- Lc
                curr.tc <- VBGF(L=curr.Lc, param = list(Linf=Linf,K=K,t0=t0,ts=ts,C=C))
                # current exploitation rate
                curr.F = (M * curr.E)/(1-curr.E)  # curr.F <- (M * curr.E)/(1-curr.E)
                tmpList <- list(Linf=Linf,
                                Winf = Winf,
                                K = K,
                                M = M,
                                t0 = t0,
                                tr = tr,
                                tc = curr.tc)

                if(length(s_list) == 1 | selecType == "knife_edge"){
                    tmpRES <- ypr(param = tmpList, FM_change = curr.F)
                }
                if(length(s_list) > 1 & selecType != "knife_edge"){
                    P <- select_ogive(s_list, Lt =  Lt, Lc = curr.Lc)
                    tmpRES <- ypr_sel(param = tmpList, FM_change = curr.F, Lt, P)
                    tmpRES$yr <- tmpRES$ryr * Winf * exp(M * (tr - t0))
                    tmpRES$br <- tmpRES$rbr * Winf * exp(M * (tr - t0))
                }

                F01[bi] <- FM_change[N01]
                Fmax[bi] <- FM_change[Nmax]

                if(length(B_R.percent) > 0){
                    F05[bi] <- FM_change[N05]
                }else{
                    F05[bi] <- NA
                }


                if("L50" %in% names(bootRaw)){
                    if(is.na(bootRaw$L50[bi]) | is.nan(bootRaw$L50[bi]) |
                       is.null(bootRaw$L50[bi]) | bootRaw$L50[bi] == Inf |
                   bootRaw$L50[bi] == -Inf | bootRaw$L50[bi] < 0 |
                   bootRaw$L75[bi] >= bootRaw$Linf[bi]){
                    F01[bi] <- NaN  ## hack
                    Fmax[bi] <- NaN
                    F05[bi] <- NaN
                    }
                }
                if(is.na(bootRaw$Z[bi]) | is.nan(bootRaw$Z[bi]) |
                   is.null(bootRaw$Z[bi]) | bootRaw$Z[bi] == Inf |
                   bootRaw$Z[bi] == -Inf){
                    F01[bi] <- NaN  ## hack
                    Fmax[bi] <- NaN
                    F05[bi] <- NaN
                }
            }
        }
        bootRaw <- cbind(bootRaw, data.frame(F01,Fmax,F05))
        tmp <- as.data.frame(bootRaw[,(ncol(boot$bootRaw)+(1:3))])


        ## estimate FF01, FFmax, FF05 if FM in dataframe
        if("FM" %in% colnames(bootRaw)){
            bootRaw[,(ncol(bootRaw)+1)] <- bootRaw$FM / F01
            bootRaw[,(ncol(bootRaw)+1)] <- bootRaw$FM / Fmax
            bootRaw[,(ncol(bootRaw)+1)] <- bootRaw$FM / F05
            colnames(bootRaw) <- c(colnames(bootRaw)[-((ncol(bootRaw)-2):ncol(bootRaw))],
                                   c("FF01","FFmax","FF05"))
            tmp <- cbind(tmp,as.data.frame(bootRaw[,((ncol(bootRaw)-2):ncol(bootRaw))]))
        }

        idx <- apply(tmp, 2, function(x) all(x == 0 | is.na(x) | x == 1))
        tmp <- tmp[,!idx]
        nx <- ncol(tmp)


        ## max density and CIS
        resMaxDen <- vector("numeric", ncol(tmp))
        ciList <- vector("list", ncol(tmp))
        for(i in seq(nx)){
            ## max densities
            x <- try(ks::kde(as.numeric(na.omit(tmp[,i]))), TRUE)
            if(class(x) != "try-error"){
                ind <- which(x$estimate > x$cont["99%"])
                resMaxDen[i] <- mean(x$eval.points[ind])
                ## confidence intervals
                CItxt <- paste0(round(5), "%")
                inCI <- rle( x$estimate > x$cont[CItxt] )
                start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
                end.idx <- cumsum(inCI$lengths)
                limCI <- try(range(x$eval.points[start.idx[min(which(inCI$values),na.rm=TRUE)]:
                                                 end.idx[max(which(inCI$values),na.rm=TRUE)]]))
                if(class(limCI) != "try-error"){  ## haven't quite figured out why limCI can give NA, but either all of x$eval.points or all of start.idx or end.idx are NA...
                    ciList[[i]] <- limCI
                }else{
                    ciList[[i]] <- c(NA,NA)
                }
            }else{
                if(all(!is.na(tmp[,i]))){
                    tmp5 <- tmp[,i]
                    tmp6 <- names(table(tmp5))[which.max(table(tmp5))]
                    resMaxDen[i] <- as.numeric(tmp6)
                    ciList[[i]] <- c(NaN,NaN)
                }else{
                    resMaxDen[i] <- NaN
                    ciList[[i]] <- c(NaN,NaN)
                }
            }
        }

        resCIs <- cbind(boot$CI,t(do.call(rbind,ciList)))
        colnames(resCIs) <- c(colnames(bootRaw)[-((ncol(bootRaw)-5):ncol(bootRaw))],
                              colnames(tmp))
        rownames(resCIs) <- c("lo","up")
        resMaxDen <- cbind(boot$maxDen, t(as.data.frame(resMaxDen)))
        colnames(resMaxDen) <- c(colnames(bootRaw)[-((ncol(bootRaw)-5):ncol(bootRaw))],
                                 colnames(tmp))
        rownames(resMaxDen) <- "maxDen"


        ret <- list()
        ret$bootRaw <- bootRaw
        ret$seed <- boot$seed
        ret$maxDen <- resMaxDen
        ret$CI <- resCIs
        ret$FMvecVPA <- boot$FMvecVPA

        class(ret) <- "lfqBoot"
        return(ret)

    }else if(!is.null(boot) & class(boot) != "lfqBoot"){
        stop("You provided an object for boot, but it does not have class 'lfqBoot'. Please check.")
    }else{

        res <- param
        res$FM_relative = FM_relative

        # Beverton and Holt's ypr
        if(type == "ypr"){
          if(FM_relative) stop(noquote("ypr does not work with relative changes in FM, please provide absolute values."))
          if("par" %in% names(res)){
              Winf <- ifelse("Winf" %in% names(res$par), res$par$Winf, NA)
              Linf <- ifelse("Linf" %in% names(res$par), res$par$Linf, NA)
              K <- res$par$K
              t0 <- ifelse("t0" %in% names(res$par), res$par$t0, 0)
              C <- ifelse("C" %in% names(res$par), res$par$C, 0)
              ts <- ifelse("ts" %in% names(res$par), res$par$ts, 0)
          }else{
              Winf <- ifelse("Winf" %in% names(res), res$Winf, NA)
              Linf <- ifelse("Linf" %in% names(res), res$Linf, NA)
              K <- res$K
              t0 <- ifelse("t0" %in% names(res), res$t0, 0)
              C <- ifelse("C" %in% names(res), res$C, 0)
              ts <- ifelse("ts" %in% names(res), res$ts, 0)
          }


          M <- res$M
          a <- res$a  # might be NULL
          b <- res$b  # might be NULL
          ## this makes sure that a is always in kg/cm3 = and therefore meanBodyWeight is always in kg
          if(LW_unit == "g"){
              a <- a / 1000
          }
          if(!is.na(Linf) & is.na(Winf) & "a" %in% names(res) & "b" %in% names(res)){
            Winf <- a * (Linf ^ b)
          }
          #Linf <- ifelse(!is.null(res$Linf),res$Linf, exp(log(Winf/a)/b))
          # REALLY ? maybe without Linf: then message that Winf has to exist
          #if(is.null(Linf) | is.na(Linf)) stop("Either Linf or Winf with a and b has to be provided!")
          #if(is.null(Winf)) Winf <-  a * Linf ^ b  ###exp((log(Linf) - a)/b) # might still be NULL
          # or              Winf <- exp(log(Linf-a)/b)

          if(length(FM_change) == 1 & is.na(FM_change[1]) &
             length(E_change) == 1 & is.na(E_change[1])){
            FM_change <- seq(0,10,0.1)
            print(noquote("No fishing mortality (FM_change) or exploitation rate (E_change) was provided, a default range for fishing mortality of 0 to 10 is used."))
          }

          # transfer E_change into F_change if provided
          # if(length(FM_change) == 1 & is.na(FM_change[1]) & length(E_change) != 1 & !is.na(E_change[1])){
          #   FM_change <- (E_change * M) / (1 - E_change)
          #   FM_change[FM_change == Inf] <- (0.9999 * M) / (1 - 0.9999)
          # }
          if(length(FM_change) == 1 & is.na(FM_change[1]) &
             length(E_change) != 1 & !is.na(E_change[1])){
            E_change <- E_change[E_change <= 0.9]
            FM_change <- (E_change * M) / (1 - E_change)
          }

          # Recruitment  - knife edge
          tr <- res$tr   # might be NULL
          Lr <- res$Lr   # might be NULL
          if(is.null(tr) & is.null(Lr)) stop("Either the age or the length at recruitment (tr or Lr) has to be provided in param!")

          if(!is.na(Linf)){
              if(is.null(tr)) tr <- VBGF(L=Lr,param = list(Linf=Linf,K=K,t0=t0,C=C,ts=ts))
              ## VBGF(L=Lr,Linf=Linf,K=K,t0=t0)
              if(is.null(Lr)) Lr <- VBGF(t=tr,param = list(Linf=Linf,K=K,t0=t0,C=C,ts=ts))
              ## VBGF(t=tr,Linf=Linf,K=K,t0=t0)
          }

          # Selectivity - knife edge or with selctivtiy ogive
          tc <- res$tc   # might be NULL
          Lc <- res$Lc   # might be NULL
          if(is.null(tc) & is.null(Lc)){
            if("L50" %in% s_list) Lc <- s_list$L50
            if("Lc" %in% s_list) Lc <- s_list$Lc
            #if(!("Lc" %in% s_list) & !("L50" %in% s_list))stop("Either the age or the length at first capture (tc or Lc) has to be provided in param! \n Or provide a Lc value in s_list!")
          }
          if(!is.na(Linf)){
            if(is.null(tc) & !is.null(Lc)) tc <- VBGF(L=Lc, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(L=Lc,Linf=Linf,K=K,t0=t0)
            if(is.null(Lc) & !is.null(tc)) Lc <- VBGF(t=tc, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(t=tc,Linf=Linf,K=K,t0=t0)
            if(is.null(tc_change) & !is.null(Lc_change)) tc_change <- VBGF(L=Lc_change, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(L=Lc_change,Linf=Linf,K=K,t0=t0)
            if(is.null(Lc_change) & !is.null(tc_change)) Lc_change <- VBGF(t=tc_change, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(t=tc_change,Linf=Linf,K=K,t0=t0)
          }
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

            # KNIFE EDGE
            if(length(s_list) == 1 ){#| selecType == "knife_edge"){
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
            if(length(s_list) > 1 ){#& selecType != "knife_edge"){
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

              input <- list(Linf = ifelse(length(Linf) > 0, Linf, NA),
                            Winf = ifelse(length(Winf) > 0, Winf, NA),
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
                                      E = E)
            if(length(Ty) > 0) results.PBH$Ty <- Ty
            if(length(Ly) > 0) results.PBH$Ly <- Ly
            if(length(Wy) > 0) results.PBH$Wy <- Wy

            results.PBH$Y_R.rel <- Y_R.rel
            results.PBH$B_R.rel <- B_R.rel

            # WHY NECESSARY???
            if(length(Y_R) > 0) results.PBH$Y_R = Y_R
            if(length(B_R) > 0) results.PBH$B_R = B_R
            if(length(B_R.percent) > 0) results.PBH$B_R.percent = B_R.percent

            list_Lc_runs[[i]] <- results.PBH


            # reference points
            Nmax <- which.max(Y_R.rel)  #  should be the same as which.min(abs(deri)) which is also labelled Nmax
            deri_pot <- deri[1:Nmax]
            N01 <- which.min(abs(deri_pot - (deri[1] * 0.1)))
            N05 <- which.min(abs(B_R.percent - 50))  #which.min(abs(deri - (deri[1] * 0.5)))

            df_loop_Es <- data.frame(Lc = ifelse(!is.null(Lci),Lci,NA),
                                     tc = ifelse(!is.null(tci),tci,NA),
                                     F01 = FM_change[N01],
                                     Fmax = FM_change[Nmax])
            if(length(B_R.percent) > 0) df_loop_Es$F05 <- FM_change[N05]   # WHY NECESSARY????
            # df_loop_Es$Fmax <- FM_change[Nmax]
            df_loop_Es$E01 <- E[N01]
            df_loop_Es$Emax <- E[Nmax]
            if(length(B_R.percent) > 0) df_loop_Es$E05 <- E[N05]    # WHY NECESSARY????
            # df_loop_Es$Emax <- E[Nmax]

            list_Es[[i]] <- df_loop_Es

            # update counter and progress bar
            if (!hide.progressbar) {
             if(nlk > 1){
              setTxtProgressBar(pb, counter)
              counter <- counter + 1
             }
            }
          }

          df_Es <- do.call(rbind,list_Es)

          names(list_Lc_runs) <- paste0("Lc_", Lc_tc)   # names(list_tc_runs) <- tc
          ret <- c(res,list(FM_change = FM_change,
                            Lc = Lc,
                            tc = tc,
                            list_Lc_runs = list_Lc_runs,   #   list_tc_runs = list_tc_runs,
                            df_Es = df_Es))   #   df_Es = df_Es,


          if(!is.na(curr.E) & !is.na(curr.Lc)){
            curr.tc <- VBGF(L=curr.Lc, param = list(Linf=Linf,K=K,t0=t0))
            # current exploitation rate
            curr.F = (M * curr.E)/(1-curr.E)  # curr.F <- (M * curr.E)/(1-curr.E)
            tmpList <- list(Linf=Linf,
                            Winf = Winf,
                            K = K,
                            M = M,
                            t0 = t0,
                            tr = tr,
                            tc = curr.tc)
            if(length(s_list) == 1 | selecType == "knife_edge"){
              tmpRES <- ypr(param = tmpList, FM_change = curr.F)
            }
            if(length(s_list) > 1 & selecType != "knife_edge"){
              P <- select_ogive(s_list, Lt =  Lt, Lc = curr.Lc)
              tmpRES <- ypr_sel(param = tmpList, FM_change = curr.F, Lt, P)
              tmpRES$yr <- tmpRES$ryr * Winf * exp(M * (tr - t0))
              tmpRES$br <- tmpRES$rbr * Winf * exp(M * (tr - t0))
            }

            df_currents <- data.frame(curr.Lc = curr.Lc,
                                      curr.tc = curr.tc,
                                      curr.E = curr.E,
                                      curr.F = curr.F,
                                      curr.YPR = tmpRES$yr,        #ypr(curr.F, curr.Lc_tc)       #, type = "length"),           # curr.YPR = ypr(curr.F, curr.Lc_tc, type = "age"),
                                      curr.YPR.rel = tmpRES$ryr,     #ypr.rel(curr.F, curr.Lc_tc),   #, type = "length"),   # curr.YPR.rel = ypr.rel(curr.F, curr.Lc_tc, type = "age"),
                                      curr.BPR = tmpRES$br,         #bpr(curr.F, curr.Lc_tc),           #, type = "length"),           # curr.BPR = bpr(curr.F, curr.Lc_tc, type = "age"),
                                      curr.BPR.rel = tmpRES$rbr)     #bpr.rel(curr.F, curr.Lc_tc))   #, type = "length"))   # curr.BPR.rel = bpr.rel(curr.F, curr.Lc_tc, type = "age"))
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

          if("par" %in% names(res)){
              Winf <- ifelse("Winf" %in% names(res$par), res$par$Winf, NA)
              Linf <- ifelse("Linf" %in% names(res$par), res$par$Linf, NA)
              K <- res$par$K
              t0 <- ifelse("t0" %in% names(res$par), res$par$t0, 0)
              C <- ifelse("C" %in% names(res$par), res$par$C, 0)
              ts <- ifelse("ts" %in% names(res$par), res$par$ts, 0)
          }else{
              Winf <- ifelse("Winf" %in% names(res), res$Winf, NA)
              Linf <- ifelse("Linf" %in% names(res), res$Linf, NA)
              K <- res$K
              t0 <- ifelse("t0" %in% names(res), res$t0, 0)
              C <- ifelse("C" %in% names(res), res$C, 0)
              ts <- ifelse("ts" %in% names(res), res$ts, 0)
          }


          # Selectivity - knife edge or with selctivtiy ogive
          tc <- res$tc   # might be NULL
          Lc <- res$Lc   # might be NULL
          if(is.null(tc) & is.null(Lc)){
            if("L50" %in% s_list) Lc <- s_list$L50
            if("Lc" %in% s_list) Lc <- s_list$Lc
            #if(!("Lc" %in% s_list) & !("L50" %in% s_list))stop("Either the age or the length at first capture (tc or Lc) has to be provided in param! \n Or provide a Lc value in s_list!")
          }
          if(!is.null(Linf)){
            if(is.null(tc) & !is.null(Lc)) tc <- VBGF(L=Lc, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(L=Lc,Linf=Linf,K=K,t0=t0)
            if(is.null(Lc) & !is.null(tc)) Lc <- VBGF(t=tc, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(t=tc,Linf=Linf,K=K,t0=t0)
            if(is.null(tc_change) & !is.null(Lc_change)) tc_change <- VBGF(L=Lc_change, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(L=Lc_change,Linf=Linf,K=K,t0=t0)
            if(is.null(Lc_change) & !is.null(tc_change)) Lc_change <- VBGF(t=tc_change, param = list(Linf=Linf,K=K,t0=t0)) # VBGF(t=tc_change,Linf=Linf,K=K,t0=t0)
          }
          tc <- c(tc,tc_change)
          Lc <- c(Lc,Lc_change)

          # age based
          if('age' %in% names(res)) classes <- as.character(res$age)
          # length based
          if('midLengths' %in% names(res)) classes <- as.character(res$midLengths)

          # create column without plus group (sign) if present
          classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
          classes.num <- as.numeric(classes.num[,1])

          if(length(FM_change) == 1 & is.na(FM_change[1]) &
             length(E_change) == 1 & is.na(E_change[1])){
            FM_change <- seq(0,10,0.1)
            print(noquote("No fishing mortality (FM_change) or exploitation rate (E_change) \nwas provided, a default range for the absolute \nfishing mortality of 0 to 10 is used."))
          }

          # transfer E_change into F_change if provided
          if(length(FM_change) == 1 & is.na(FM_change[1]) &
             length(E_change) != 1 & !is.na(E_change[1])){
            E_change <- E_change[E_change <= 0.9]
            FM_change <- (E_change * nM) / (1 - E_change)
          }
          if(length(E_change) == 1 & is.na(E_change[1])){
            E_change <- FM_change / (FM_change + nM)
          }

          Lt <- classes.num  # if age in names(res) Lt here is age because classes.num = age


          # Only FM change provided without Lc_tc change
          if((is.null(tc_change) & is.null(Lc_change))){  #  | length(s_list) == 1){

            #if(is.null(res$FM) | length(res$FM) == 1) stop(noquote("Please provide fishing mortality FM (in 'param') as a vector per size class!"))

            if(is.null(res$FM)) stop(noquote("Please provide fishing mortality FM (in 'param')!"))
            if(length(res$FM) == 1){
              if(length(s_list) > 1 | !is.null(Lc[1])){
                ##print(noquote("Fishing mortality per length class not povided, using selectivity information to derive fishing mortality per length class."))
                if(length(s_list) == 1){
                  s_list <- list(selecType = "knife_edge", L50 = Lc[1])
                }
                sel <- select_ogive(s_list, Lt = Lt)
                FM <- res$FM * sel
              }else{
                stop(noquote("Please provide either fishing mortality FM (in 'param') per length class or a Lc value!"))
              }
            }

            #prediction based on f_change
            if(!FM_relative){
              pred_mat <- as.matrix(FM/max(FM, na.rm = TRUE)) %*% FM_change
            }
            if(FM_relative){
              pred_mat <- as.matrix(FM) %*% FM_change
            }


            pred_res_list <- list()
            for(x7 in 1:length(FM_change)){
              param$Z <- pred_mat[,x7] + nM
              param$FM <- pred_mat[,x7]
              resL <- stock_sim(param, age_unit = age_unit,
                                stock_size_1 = stock_size_1, plus_group = plus_group)
              pred_res_list[[x7]] <- resL$totals
            }


            pred_res_df <- do.call(rbind, pred_res_list)
            pred_res_df$FM_change <- FM_change
            pred_res_df$E_change <- E_change

            res2 <- pred_res_df
            res3 <- c(res,res2)

              
              ## reference points
              ## Fmax
              Nmax <- which.max(pred_res_df$totY)
              ## F01 (proxy for Fmsy in ICES)
              slopeOrg <- (pred_res_df$totY[2] - pred_res_df$totY[1]) / (FM_change[2] - FM_change[1])
              slope01 <- round(0.1*slopeOrg, 2)
              slopes <- rep(NA, length(FM_change))
              slopes[1] <- slopeOrg
              for(i in 3:length(FM_change)){
                  slopes[i-1] <- round((pred_res_df$totY[i] - pred_res_df$totY[i-1]) /
                      (FM_change[i] - FM_change[i-1]),2)
              }
              dif <- abs(slopes - slope01)
              dif[is.na(dif)] <- 1e+11
              difpot <- dif[1:Nmax]
              N01 <- which.min(difpot)
              ## F05
              Bper <- rep(NA,length(pred_res_df$meanB))
              Bper[1] <- 100
              for(ix in 2:length(Bper)){
                Bper[ix] <- pred_res_df$meanB[ix]/pred_res_df$meanB[1] * 100
              }
              N05 <- which.min(abs(Bper - 50))

              if(!is.null(Lc[1]) & !is.null(tc[1])){
                df_Es <- data.frame(Lc = Lc,
                                    tc = tc,
                                    F01 = FM_change[N01],                                    
                                    Fmax = FM_change[Nmax],
                                    F05 = FM_change[N05],
                                    E01 = E_change[N01],                                    
                                    Emax = E_change[Nmax],
                                    E05 = E_change[N05])
              }else{
                  df_Es <- data.frame(F01 = FM_change[N01],
                      Fmax = FM_change[Nmax],
                      F05 = FM_change[N05],
                                    E01 = E_change[N01],                      
                                    Emax = E_change[Nmax],
                                    E05 = E_change[N05])
              }


              ret <- c(res3, list(df_Es = df_Es))


            if(!is.na(curr.E)){
              if(!is.na(curr.Lc)){
                curr.tc <- VBGF(L=curr.Lc, param = list(Linf=Linf, K=K, t0=t0))
              }else curr.tc <- NA
              # current exploitation rate
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
              res2 <- stock_sim(param = param.loop, age_unit = age_unit,
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

          # FM and Lc_tc change provided
          if(!is.null(tc_change) | !is.null(Lc_change)){
            # instead of s_list the outcome of one of the other select functions?

            if(length(s_list) == 1){
              s_list <- list(selecType = "knife_edge", L50 = Lc[1])
            }
            sel <- select_ogive(s_list, Lt = Lt) #classes.num

            sel.list <- list()
            for(x19 in 1:length(Lc)){
              sel.list[[x19]] <- select_ogive(s_list, Lt = Lt, Lc = Lc[x19]) #classes.num
            }
            Lc_mat <- do.call(cbind,sel.list)
            colnames(Lc_mat) <- Lc

            Lc_mat_FM <- Lc_mat   #max(FM, na.rm=TRUE)  # with one it should correspond to actual fishing mortality not to change in mortality (x factor)

            #list with FM_Lc_matrices per FM_change
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

            for(x21 in 1:length(FM_Lc_com_mat.list)){  #loop for length of list == Lc changes
              mati <- FM_Lc_com_mat.list[[x21]]

              pred.FM_Lc_com_res_loop1_list <- list()
              for(x22 in 1:dim(mati)[2]){

                param.loop$FM <- mati[,x22]
                param.loop$Z <- mati[,x22] + nM
                res2 <- stock_sim(param = param.loop, age_unit = age_unit,
                                  stock_size_1 = stock_size_1, plus_group=plus_group)
                pred.FM_Lc_com_res_loop1_list[[x22]] <- res2$totals

                # update counter and progress bar
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

            #for catch
            mat_FM_Lc_com.C <- do.call(rbind, pred.FM_Lc_com_res_loopC_list)
            rownames(mat_FM_Lc_com.C) <- Lc
            colnames(mat_FM_Lc_com.C) <- FM_change

            #for yield
            mat_FM_Lc_com.Y <- do.call(rbind, pred.FM_Lc_com_res_loopY_list)
            rownames(mat_FM_Lc_com.Y) <- Lc
            colnames(mat_FM_Lc_com.Y) <- FM_change

            #for biomass
            mat_FM_Lc_com.B <- do.call(rbind, pred.FM_Lc_com_res_loopB_list)
            rownames(mat_FM_Lc_com.B) <- Lc
            colnames(mat_FM_Lc_com.B) <- FM_change

            #for value
            mat_FM_Lc_com.V <- do.call(rbind, pred.FM_Lc_com_res_loopV_list)
            rownames(mat_FM_Lc_com.V) <- Lc
            colnames(mat_FM_Lc_com.V) <- FM_change

            # transvers matrices for plotting (the opposite arrangement from book)
            mat_FM_Lc_com.C <- t(mat_FM_Lc_com.C)
            mat_FM_Lc_com.Y <- t(mat_FM_Lc_com.Y)
            mat_FM_Lc_com.B <- t(mat_FM_Lc_com.B)
            mat_FM_Lc_com.V <- t(mat_FM_Lc_com.V)


              ## reference points
              ## F01 (proxy for Fmsy in ICES)
              slopes <- matrix(NA,ncol=dim(mat_FM_Lc_com.Y)[2],
                               nrow=dim(mat_FM_Lc_com.Y)[1])
              slopeOrg <- (mat_FM_Lc_com.Y[2,] - mat_FM_Lc_com.Y[1,]) /
                  (rep(FM_change[2],dim(mat_FM_Lc_com.Y)[2]) -
                   rep(FM_change[1],dim(mat_FM_Lc_com.Y)[2]))
              slope01 <- round(0.1*slopeOrg, 2)
              slopes[1,] <- slopeOrg
              for(i in 3:length(FM_change)){
                  slopes[i-1,] <- round((mat_FM_Lc_com.Y[i,] - mat_FM_Lc_com.Y[i-1,]) /
                                        (rep(FM_change[i],dim(mat_FM_Lc_com.Y)[2]) -
                                         rep(FM_change[i-1],dim(mat_FM_Lc_com.Y)[2])),2)
              }
              dif <- t(apply(slopes,1,function(x) abs(x - slope01)))
              dif[is.na(dif)] <- 1e+11
              difpot <- dif[1:Nmax,]                                          ### DOUBLE CHECK
              N01 <- apply(difpot, MARGIN = 2, FUN = which.min)
              ## F05
              mat_FM_Lc_com.Bper <- matrix(NA,ncol=dim(mat_FM_Lc_com.B)[2],
                                                   nrow=dim(mat_FM_Lc_com.B)[1])
              mat_FM_Lc_com.Bper[1,] <- 100
              for(ix in 2:dim(mat_FM_Lc_com.B)[1]){
                mat_FM_Lc_com.Bper[ix,] <- mat_FM_Lc_com.B[ix,]/mat_FM_Lc_com.B[1,] *100
              }
              N05 <- apply(mat_FM_Lc_com.Bper, MARGIN = 2,
                           FUN = function(x) which.min(abs(x - 50)))
              ## Fmax
              Nmax <- apply(mat_FM_Lc_com.Y, MARGIN = 2, FUN = which.max)

              if((!is.null(Lc[1]) & !is.null(tc[1])) | (!is.na(Lc[1]) & !is.na(tc[1])) ){
                df_Es <- data.frame(Lc = Lc,
                                    tc = tc,
                                    F01 = FM_change[N01],
                                    Fmax = FM_change[Nmax],
                                    F05 = FM_change[N05],
                                    E01 = E_change[N01],                                    
                                    Emax = E_change[Nmax],
                                    E05 = E_change[N05])
              }else{
                  df_Es <- data.frame(
                      F01 = FM_change[N01],                      
                      Fmax = FM_change[Nmax],
                      F05 = FM_change[N05],
                      E01 = E_change[N01],
                      Emax = E_change[Nmax],
                      E05 = E_change[N05])
              }


            ret <- c(res,
                     list(FM_change = FM_change,
                          # FM_relative = FM_relative,
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
                curr.tc <- VBGF(L=curr.Lc, param = list(Linf=Linf, K=K, t0=t0))
              }else curr.tc <- NA

              # current exploitation rate
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

        # return results and plot
        class(ret) <- "predict_mod"
        if(plot) plot(ret, mark = mark)

        return(ret)
    }

}
