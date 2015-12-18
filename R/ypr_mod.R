#' @title Prediction models: Beverton and Holt's yield per recruit model
#
#' @description  This is a function ..
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{Winf} or \strong{Linf}: infinite weight (or length) for investigated
#'   species in cm [cm] or gramm [g],
#'   \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'   \item \strong{t0}: theoretical time zero, at which individuals of this species hatch (default: 0),
#'   \item \strong{M}: natural mortality [1/year],
#'   \item \strong{a}: length-weight relationship coefficent (W = a * L^b),
#'   \item \strong{b}: length-weight relationship coefficent (W = a * L^b),
#'   \item \strong{tr} or \strong{Lr}: age or length of recruitment to fishery, respectively;
#' }
#' @param tc_Lc age or length at first capture
#' @param FM  fishing mortality (Default: FM = NA)
#'
#' @keywords function prediction yield-per-recruit
#'
#' @examples
#' #______________________________________
#' # age structured data
#' # Nemipterus marginatus
#' threadfin <- list(Winf = 286,K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4)
#'
#' # run model
#' ypr_mod(threadfin, tc_Lc = seq(0.2,1,0.2), FM = seq(0,6,0.1))  #where it is maximal  = MSY
#'
#' # Leiognathus spendens (Pauly 1980)
#' ponyfish <- list(Winf = 64, K = 1, t0 = -0.2, M = 1.8, tr = 0.2)
#'
#' # run model
#' ypr_mod(ponyfish, tc_Lc = c(0.2,0.3,1.0))
#'
#' #______________________________________
#' # length structured data
#' # Xiphias gladius (Berkeley and Houde 1980)
#' swordfish <- list(Linf = 309, K = 0.0949, M = 0.18, a=0.0003, b=3, Lr = 90)  ## T_Lr , a, b ??? assumed
#'
#' # run model
#' ypr_mod(tc_Lc = c(100,118,150,180))
#'
#'
#'
#'
#' ####test: E <- seq(0,0.9,0.1) FM <- E * M / (1 -E)
#'
#' @details ...
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{tplusdt_2} or \strong{t_midL}: relative,
#'   \item \strong{lnC_dt}: rearranged,
#'   \item \strong{reg_int}: the,
#'   \item \strong{Z}: the,
#'   \item \strong{se}: the,
#'   \item \strong{intercept}: intercep;
#' }
#'
#' @references
#' Berkeley and Houde 1980.
#'
#' Pauly, D., 1980.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

param = threadfin


ypr_mod <- function(param, tc_Lc, FM = NA){

  res <- param
  M <- res$M
  K <- res$K
  t0 <- ifelse(!is.null(res$t0),res$t0,0)

  if(length(FM) == 1 & is.na(FM[1])){
    FM <- seq(0,10,0.1)
    warning("No fishing mortality (FM) was provided, a default range of 0 to 10 is used.")
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                        Age data                          #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if("Winf" %in% names(res)){
    Winf <- res$Winf
    tr <- res$tr
    tc <- tc_Lc
    list_tc_runs <- list()
    for(i in 1:length(tc)){
      tci <- tc[i]

      #Biomass per Recruit
      S <- exp(-K * (tci - t0))
      B_R <- exp(-M*(tci-tr)) * Winf *
        ((1/(FM+M)) - ((3*S)/((M+FM)+K)) +
           ((3*(S^2))/((M+FM)+(2*K))) - ((S^3)/((M+FM) + (3*K))))

      #virgin biomass
      Bv_R <- B_R[which(FM == 0)]
      #biomass of exploited part of the cohort (biomass of fish older than tc)

      #biomass in percetage of virgin biomass
      B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

      #Yield per Recruit
      Y_R <- B_R * FM

      #relative yield per recruit - mostly done with length frequency data (exclusively?)
      Y_R.rel <- Y_R * (exp(-M *(tr - t0))) / Winf


      #mean age in annual yield
      Ty <- (1 / (M+FM)) + tci

      #mean length in the annual yield
      #Ly <- Linf * (1 - (((M+FM)*S)/((M+FM)+K)))

      #mean weight in annual yield
      Wy <- (M+FM) * Winf *
        ((1/(FM+M)) - ((3*S)/((M+FM)+K)) +
           ((3*(S^2))/((M+FM)+(2*K))) - ((S^3)/((M+FM) + (3*K))))


      results.PBH <- data.frame(FM = FM,
                                Y_R = Y_R,
                                Y_R.rel = Y_R.rel,
                                B_R = B_R,
                                B_R.percent = B_R.percent,
                                Ty = Ty,
                                Wy = Wy)


      list_tc_runs[[i]] <- results.PBH

    }
    names(list_tc_runs) <- tc
    ret <- list_tc_runs
    class(ret) <- "ypr_mod"
    # plot results
    plot(ret)
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if("Linf" %in% names(res)){
    Linf <- res$Linf
    Lr <- res$Lr
    Lc <- tc_Lc
    list_Lc_runs <- list()
    for(i in 1:length(Lc)){

      Lci <- Lc[i]

      E <- FM/(FM + M)
      # transform Linf in Winf #### CHECK!!
      Winf <- a * (Linf ^ b)

      #yiel per recurit for length data   # also possbile inputing option: F/K
      S <- 1 - (Lci/Linf)  # == U  ##(U <- 1 - (Lci/Linf))
      A.PBH <- ((Linf - Lci)/(Linf-Lr))^(M/K)
      Y_R <- FM * A.PBH * Winf * ((1/(M+FM)) - (3*S)/((M+FM)+K) +
                                           (3*S^2)/((M+FM)+2*K) -
                                           (S^3)/((M+FM)+3*K))

      #biomass per recruit for length data?
      B_R <- Y_R / FM

      #virgin biomass
      Bv_R <- B_R[which(FM == 0)]
      #biomass of exploited part of the cohort (biomass of fish older than tc)

      #biomass in percetage of virgin biomass
      B_R.percent <- round((B_R / Bv_R ) * 100, digits = 1)

      #relative yield per recruit - mostly done with length frequency data (exclusively?)
      m <- K/(M+FM)
      Y_R.rel <- E * S^(M/K) * (1 - ((3*S)/(1+m)) +
                                          ((3*S^2)/(1+2*m)) - ((S^3)/(1+3*m)))

      #mean length in the annual yield
      Ly <- Linf * (1 - (((M+FM)*S)/((M+FM)+K)))

      #mean weight in annual yield
      # Wy <- (M+FM) * Winf *
      #    ((1/(FM+M)) - ((3*S)/((M+FM)+K)) +
      #       ((3*(S^2))/((M+FM)+(2*K))) - ((S^3)/((M+FM) + (3*K))))


      results.PBH <- data.frame(FM = FM,
                                Ly = Ly,
                                E = E,
                                Y_R = Y_R,
                                Y_R.rel = Y_R.PBH.rel,
                                B_R = B_R,
                                B_R.pecent = B_R.percent)


      list_Lc_runs[[i]] <- results.PBH

    }
    names(list_Lc_runs) <- Lc
    ret <- list_Lc_runs
    class(ret) <- "ypr_mod"

    # plot results
    plot(ret)
  }
}

## problem of two cases: tc and Co are given or Lc and Co. case dependent or
##different functions?
