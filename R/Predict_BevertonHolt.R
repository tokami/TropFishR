#' @title Beverton and Holt's prediction model
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param classes Midpoints of the length class as vector (length frequency data) or ages as vector (age composition data).
#' @param catch Catch as vector, or a matrix with catches of subsequent years if the catch curve with constat time intervals should be applied.
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch (default: 0).
#' @param Winf Infinite weight in grams!
#'
#' @examples
#'
#'
#' @details For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references
#'
#'
#' @export


#age data example 1
# Nemipterus marginatus
K = 0.37
M = 1.1
Tc = c(0.2,0.3,1.0)
Tr = 0.4
t0 = -0.2
Winf = 286
FM = seq(0,6,0.1)

Predict_BevertonHolt(K = 0.37,
                     M = 1.1,
                     Tc = 1.0,
                     Tr = 0.4,
                     t0 = -0.2,
                     Winf = 286,
                     FM = seq(0,6,0.1))
#where it is maximal  = MSY

############

#age data example 2
#Leiognathus spendens (Pauly 1980)
Winf = 64
K = 1
t0 = -0.2
Tr = 0.2
M = 1.8
Tc = c(0.2,0.3,1.0)

############

#length data example
#Xiphias gladius (Berkeley and Houde 1980)
Linf = 309
K = 0.0949
M = 0.18
Lc = c(100,118,150,180) #Lc = c(118,150)
Lr = 90 #??? assumed


Predict_BevertonHolt <- function(Winf = NA, Linf = NA, K, t0 = NA,
                                  M, Tr = NA, Tc = NA, Lc = NA, FM = NA){

  if(is.na(Winf) & is.na(Linf)) stop("You have to provide Linf or Winf!")


  F_PBH <- FM
  if(length(F_PBH) == 1 & is.na(F_PBH[1])){
    F_PBH <- seq(0,10,0.1)
    warning("No fishing mortality (FM) was provided, a default range of 0 to 10 is used.")
    }


  #with age data
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

    #mean age in annual yield
    Ty <- (1 / (M+F_PBH)) + Tci

    #mean length in the annual yield
    Ly <- Linf * (1 - (((M+F_PBH)*S.PBH)/((M+F_PBH)+K)))

    #mean weight in annual yield
    Wy <- (M+F_PBH) * Winf *
      ((1/(F_PBH+M)) - ((3*S.PBH)/((M+F_PBH)+K)) +
         ((3*(S.PBH^2))/((M+F_PBH)+(2*K))) - ((S.PBH^3)/((M+F_PBH) + (3*K))))


    results.PBH <- data.frame(FM = F_PBH,
                              Y_R = Y_R.PBH,
                              B_R = B_R.PBH,
                              B_R.percent = B_R.PBH.percent,
                              Ty = Ty,
                              Ly = Ly,
                              Wy = Wy)


    list_Tc_runs[[i]] <- results.PBH

  }
  names(list_Tc_runs) <- Tc

  #plot
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


  #with length data
  list_Lc_runs <- list()
  for(i in 1:length(Lc)){

    Lci <- Lc[i]

    #test:
    E <- seq(0,0.9,0.1)
    F_PBH <- E * M / (1 -E)

    #E <- F_PBH/(F_PBH + M)


    #biomass per recruit for length data?

    #yiel per recurit for length data   # also possbile inputing option: F/K
    S.PBH <- 1 - (Lci/Linf)  # == U  ##(U <- 1 - (Lci/Linf))
    A.PBH <- ((Linf - Lci)/(Linf-Lr))^(M/K)
    Y_R.PBH <- F_PBH * A.PBH * Winf * ( (1/(M+F_PBH)) - (3*S.PBH)/((M+F_PBH)+K) +
                                          (3*S.PBH^2)/((M+F_PBH)+2*K) -
                                          (S.PBH^3)/((M+F_PBH)+3*K))

    #relative yield per recruit - mostly done with length frequency data (exclusively?)
    m <- K/(M+F_PBH)
    Y_R.PBH.rel <- E * S.PBH^(M/K) * (1 - ((3*S.PBH)/(1+m)) +
                                    ((3*S.PBH^2)/(1+2*m)) - ((S.PBH^3)/(1+3*m)))

    #mean length in the annual yield
    Ly <- Linf * (1 - (((M+F_PBH)*S.PBH)/((M+F_PBH)+K)))

    #mean weight in annual yield
    Wy <- (M+F_PBH) * Winf *
      ((1/(F_PBH+M)) - ((3*S.PBH)/((M+F_PBH)+K)) +
         ((3*(S.PBH^2))/((M+F_PBH)+(2*K))) - ((S.PBH^3)/((M+F_PBH) + (3*K))))


    results.PBH <- data.frame(FM = F_PBH,
                              Ly = Ly,
                              Wy = Wy,
                              E = E,
                              Y_R.rel = Y_R.PBH.rel)


    list_Lc_runs[[i]] <- results.PBH

  }
  names(list_Lc_runs) <- Lc


  #plot
  Lc_start <- which(Lc == max(Lc,na.rm=T))
  offset_text <- list_Lc_runs[[Lc_start]]$Y_R.rel[length(list_Lc_runs[[Lc_start]]$Y_R.rel)] *
    0.02
  par(mar = c(5, 4, 4, 4) + 0.3)
  plot(list_Lc_runs[[Lc_start]]$FM, list_Lc_runs[[Lc_start]]$Y_R, type = 'l',
       ylab = "rel. Y/R", xlab = "F",lty=Lc_start)
  text(x = list_Lc_runs[[Lc_start]]$FM[length(list_Lc_runs[[Lc_start]]$FM)],
       y = (list_Lc_runs[[Lc_start]]$Y_R[length(list_Lc_runs[[Lc_start]]$Y_R)] +
              offset_text),
       labels = bquote(Lc[.(names(list_Lc_runs)[Lc_start])]))
  seq_Lc <- 1:length(list_Lc_runs)
  seq_Lc <- seq_Lc[-Lc_start]
  for(j in seq_Lc){
    lines(list_Lc_runs[[j]]$FM, list_Lc_runs[[j]]$Y_R, type = 'l',lty = j)
    text(x = list_Lc_runs[[j]]$FM[length(list_Lc_runs[[j]]$FM)],
         y = (list_Lc_runs[[j]]$Y_R[length(list_Lc_runs[[j]]$Y_R)] +
                offset_text),
         labels = bquote(Lc[.(names(list_Lc_runs)[j])]))
  } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area

} ## problem of two cases: Tc and Co are given or Lc and Co. case dependent or different functions?




