#' @title Beverton and Holt's prediction model
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param W_Linf Infinite length (or weight) for investigated species in cm [cm] or gramm [g].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch (default: 0).
#' @param M Natural mortality
#' @param T_Lr Age or Length of recruitment
#' @param T_Lc Age or length of first capture
#' @param FM  fishing mortlality Default FM = NA
#' @param a length weight relationship coefficient 1 Default = NA
#' @param b length weight relationship coefficient 2 Default = NA

#' @examples
#' # age structured data
#' # example 1 - Nemipterus marginatus
#' Predict_BevertonHolt(K = 0.37,M = 1.1,T_Lc = seq(0.2,1,0.2),
#'    T_Lr = 0.4,t0 = -0.2,W_Linf = 286,FM = seq(0,6,0.1),datatype = 'age')
#' #where it is maximal  = MSY
#'
#'
#' # example 2 - Leiognathus spendens (Pauly 1980)
#' Predict_BevertonHolt(W_Linf = 64,K = 1,t0 = -0.2,T_Lr = 0.2,M = 1.8,
#'    T_Lc = c(0.2,0.3,1.0),datatype = 'age')
#'
#' # length structured data
#' # Xiphias gladius (Berkeley and Houde 1980)
#' Predict_BevertonHolt(W_Linf = 309,K = 0.0949,M = 0.18,
#'    T_Lc = c(100,118,150,180),T_Lr = 90 ,datatype = 'length',a=0.0003,b=3)
#' ## T_Lr , a, b ??? assumed
#'
#' ####test: E <- seq(0,0.9,0.1) F_PBH <- E * M / (1 -E)
#'
#' @details For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references xxx
#'
#'
#' @export

Predict_BevertonHolt <- function(W_Linf, K, M, T_Lr, T_Lc, t0 = NA,
                                  FM = NA, datatype, a=NA, b=NA){

  F_PBH <- FM
  if(length(F_PBH) == 1 & is.na(F_PBH[1])){
    F_PBH <- seq(0,10,0.1)
    warning("No fishing mortality (FM) was provided, a default range of 0 to 10 is used.")
    }


  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                        Age data                          #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(datatype == 'age'){
    Winf <- W_Linf
    Tr <- T_Lr
    Tc <- T_Lc
    list_Tc_runs <- list()
    for(i in 1:length(Tc)){
      Tci <- Tc[i]

      #Biomass per Recruit
      S.PBH <- exp(-K * (Tci - t0))
      B_R.PBH <- exp(-M*(Tci-Tr)) * Winf *
        ((1/(F_PBH+M)) - ((3*S.PBH)/((M+F_PBH)+K)) +
           ((3*(S.PBH^2))/((M+F_PBH)+(2*K))) - ((S.PBH^3)/((M+F_PBH) + (3*K))))

      #virgin biomass
      Bv_R.PBH <- B_R.PBH[which(F_PBH == 0)]
      #biomass of exploited part of the cohort (biomass of fish older than Tc)

      #biomass in percetage of virgin biomass
      B_R.PBH.percent <- round((B_R.PBH / Bv_R.PBH ) * 100, digits = 1)

      #Yield per Recruit
      Y_R.PBH <- B_R.PBH * F_PBH

      #relative yield per recruit - mostly done with length frequency data (exclusively?)
      Y_R.PBH.rel <- Y_R.PBH * (exp(-M *(Tr - t0))) / Winf


      #mean age in annual yield
      Ty <- (1 / (M+F_PBH)) + Tci

      #mean length in the annual yield
      #Ly <- Linf * (1 - (((M+F_PBH)*S.PBH)/((M+F_PBH)+K)))

      #mean weight in annual yield
      Wy <- (M+F_PBH) * Winf *
        ((1/(F_PBH+M)) - ((3*S.PBH)/((M+F_PBH)+K)) +
           ((3*(S.PBH^2))/((M+F_PBH)+(2*K))) - ((S.PBH^3)/((M+F_PBH) + (3*K))))


      results.PBH <- data.frame(FM = F_PBH,
                                Y_R = Y_R.PBH,
                                Y_R.rel = Y_R.PBH.rel,
                                B_R = B_R.PBH,
                                B_R.percent = B_R.PBH.percent,
                                Ty = Ty,
                                Wy = Wy)


      list_Tc_runs[[i]] <- results.PBH

    }
    names(list_Tc_runs) <- Tc

    #plot Y/R & B/R
    Tc_start <- which(Tc == max(Tc,na.rm=T))
    offset_text <- list_Tc_runs[[Tc_start]]$Y_R[length(list_Tc_runs[[Tc_start]]$Y_R)] *
      0.02
    par(mar = c(5, 4, 4, 4) + 0.3)
    plot(list_Tc_runs[[Tc_start]]$FM, list_Tc_runs[[Tc_start]]$Y_R, type = 'l',
         ylab = "Y/R", xlab = "F",lty=Tc_start)
    text(x = list_Tc_runs[[Tc_start]]$FM[length(list_Tc_runs[[Tc_start]]$FM)],
         y = (list_Tc_runs[[Tc_start]]$Y_R[length(list_Tc_runs[[Tc_start]]$Y_R)] +
                offset_text),
         labels = bquote(Tc[.(names(list_Tc_runs)[Tc_start])]))
    seq_Tc <- 1:length(list_Tc_runs)
    seq_Tc <- seq_Tc[-Tc_start]
    for(j in seq_Tc){
      lines(list_Tc_runs[[j]]$FM, list_Tc_runs[[j]]$Y_R, type = 'l',
            ylab = "Y/R", xlab = "F",lty = j)
      text(x = list_Tc_runs[[j]]$FM[length(list_Tc_runs[[j]]$FM)],
           y = (list_Tc_runs[[j]]$Y_R[length(list_Tc_runs[[j]]$Y_R)] +
                  offset_text),
           labels = bquote(Tc[.(names(list_Tc_runs)[j])]))
    } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area
    par(new=T)
    plot(list_Tc_runs[[1]]$FM, list_Tc_runs[[1]]$B_R, type = 'l',
         axes = F, ylab = '', xlab ='', lty = Tc_start,
         col = 'blue')
    axis(side = 4, at = pretty(range(list_Tc_runs[[1]]$B_R)), col = 'blue',
         col.axis = 'blue')
    mtext(side = 4, text = "B/R", line = 3, col = 'blue')
    for(j in seq_Tc){
      lines(list_Tc_runs[[j]]$FM, list_Tc_runs[[j]]$B_R, type = 'l',
            ylab = "Y/R", xlab = "F", col='blue', lty = j)
    } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area
  }


  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #                       Length data                        #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(datatype == 'length'){
    Linf <- W_Linf
    Lr <- T_Lr
    Lc <- T_Lc
    list_Lc_runs <- list()
    for(i in 1:length(Lc)){

      Lci <- Lc[i]

      E <- F_PBH/(F_PBH + M)
      # transform Linf in Winf #### CHECK!!
      Winf <- a * (Linf ^ b)

      #yiel per recurit for length data   # also possbile inputing option: F/K
      S.PBH <- 1 - (Lci/Linf)  # == U  ##(U <- 1 - (Lci/Linf))
      A.PBH <- ((Linf - Lci)/(Linf-Lr))^(M/K)
      Y_R.PBH <- F_PBH * A.PBH * Winf * ((1/(M+F_PBH)) - (3*S.PBH)/((M+F_PBH)+K) +
                                            (3*S.PBH^2)/((M+F_PBH)+2*K) -
                                            (S.PBH^3)/((M+F_PBH)+3*K))

      #biomass per recruit for length data?
      B_R.PBH <- Y_R.PBH / F_PBH

      #virgin biomass
      Bv_R.PBH <- B_R.PBH[which(F_PBH == 0)]
      #biomass of exploited part of the cohort (biomass of fish older than Tc)

      #biomass in percetage of virgin biomass
      B_R.PBH.percent <- round((B_R.PBH / Bv_R.PBH ) * 100, digits = 1)

      #relative yield per recruit - mostly done with length frequency data (exclusively?)
      m <- K/(M+F_PBH)
      Y_R.PBH.rel <- E * S.PBH^(M/K) * (1 - ((3*S.PBH)/(1+m)) +
                                          ((3*S.PBH^2)/(1+2*m)) - ((S.PBH^3)/(1+3*m)))

      #mean length in the annual yield
      Ly <- Linf * (1 - (((M+F_PBH)*S.PBH)/((M+F_PBH)+K)))

      #mean weight in annual yield
     # Wy <- (M+F_PBH) * Winf *
     #    ((1/(F_PBH+M)) - ((3*S.PBH)/((M+F_PBH)+K)) +
     #       ((3*(S.PBH^2))/((M+F_PBH)+(2*K))) - ((S.PBH^3)/((M+F_PBH) + (3*K))))


      results.PBH <- data.frame(FM = F_PBH,
                                Ly = Ly,
                                E = E,
                                Y_R = Y_R.PBH,
                                Y_R.rel = Y_R.PBH.rel,
                                B_R = B_R.PBH,
                                B_R.pecent = B_R.PBH.percent)


      list_Lc_runs[[i]] <- results.PBH

    }
    names(list_Lc_runs) <- Lc

    #plot Y/R & B/R
    Lc_start <- which(Lc == max(Lc,na.rm=T))
    offset_text <- list_Lc_runs[[Lc_start]]$Y_R[length(list_Lc_runs[[Lc_start]]$Y_R)] *
      0.02
    par(mar = c(5, 4, 4, 4) + 0.3)
    plot(list_Lc_runs[[Lc_start]]$FM, list_Lc_runs[[Lc_start]]$Y_R, type = 'l',
         ylab = "Y/R", xlab = "F",lty=Lc_start)
    text(x = list_Lc_runs[[Lc_start]]$FM[length(list_Lc_runs[[Lc_start]]$FM)],
         y = (list_Lc_runs[[Lc_start]]$Y_R[length(list_Lc_runs[[Lc_start]]$Y_R)] +
                offset_text),
         labels = bquote(Lc[.(names(list_Lc_runs)[Lc_start])]))
    seq_Lc <- 1:length(list_Lc_runs)
    seq_Lc <- seq_Lc[-Lc_start]
    for(j in seq_Lc){
      lines(list_Lc_runs[[j]]$FM, list_Lc_runs[[j]]$Y_R, type = 'l',
            ylab = "Y/R", xlab = "F",lty = j)
      text(x = list_Lc_runs[[j]]$FM[length(list_Lc_runs[[j]]$FM)],
           y = (list_Lc_runs[[j]]$Y_R[length(list_Lc_runs[[j]]$Y_R)] +
                  offset_text),
           labels = bquote(Lc[.(names(list_Lc_runs)[j])]))
    } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area
    par(new=T)
    plot(list_Lc_runs[[1]]$FM, list_Lc_runs[[1]]$B_R, type = 'l',
         axes = F, ylab = '', xlab ='', lty = Lc_start,
         col = 'blue')
    axis(side = 4, at = pretty(range(list_Lc_runs[[1]]$B_R)), col = 'blue',
         col.axis = 'blue')
    mtext(side = 4, text = "B/R", line = 3, col = 'blue')
    for(j in seq_Lc){
      lines(list_Lc_runs[[j]]$FM, list_Lc_runs[[j]]$B_R, type = 'l',
            ylab = "Y/R", xlab = "F", col='blue', lty = j)
    } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area
  }
} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?
