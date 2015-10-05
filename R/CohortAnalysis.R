#' @title Cohort analysis
#
#' @description Cohort analysis
#'
#' @param midLengths Midpoints of the length class as vector
#'
#' @param catch Catch per sampling time as matrix or the total catch as vector.
#'
#' @param Linf Infinite length for investigated species in cm [cm].
#'
#' @param K Growth coefficient for investigated species per year [1/year].
#'
#' @param t0 Theoretical time zero, at which individuals of this species hatch.
#'
#' @param catchCorFac optional: Correction factor for catch, in case provided catch does spatially or temporarily not reflect catch for fishing ground of a whole year.
#'
#' @param M Natural mortality [1/year]
#'
#' @param terminalF terminal fishing mortality
#'
#' @param a length-weight relationship coefficent (W = a * L^b)
#'
#' @param b length-weight relationship coefficent (W = a * L^b)
#'
#' @examples
#' data("ex.CohortAnalysis")
#' output = CohortAnalysis(midLengths = ex.CohortAnalysis$midLengths, catch = ex.CohortAnalysis$catch, Linf = 130, K = 0.1, M = 0.28, terminalF = 0.28, a = 0.00001, b = 3)
#' output
#' @details Cohort analysis
#'
#' @export

CohortAnalysis <- function(midLengths, catch, Linf, K, t0 = 0, M, terminalF,
                           a, b, catchCorFac = NA){

  if(length(midLengths) != length(catch)) stop("midLengths and catch do not have the same length!")

  df.CA <- data.frame(midLengths.CA = midLengths,
                      catch.CA = catch)

  #calculate size class interval
  interval.CA <- df.CA$midLengths.CA[2] - df.CA$midLengths.CA[1]

  # t of lower length classes
  df.CA$lowerLength.CA <- df.CA$midLengths.CA - (interval.CA / 2)
  df.CA$catchCor.CA <- df.CA$catch.CA * catchCorFac
  if(is.na(catchCorFac)) df.CA$catchCor.CA <- df.CA$catch.CA
  df.CA$t_L1 <- (t0 - (1/K)) * log(1 - (df.CA$lowerLength.CA / Linf))

  # delta t
  df.CA$dt <- NA
  for(x1 in 1:(length(df.CA$dt)-1)){
    df.CA$dt[x1] <- df.CA$t_L1[x1+1] - df.CA$t_L1[x1]
  }

  # t of midlengths
  df.CA$t_midL <- (t0 - (1/K)) * log(1 - (df.CA$midLengths.CA / Linf))

  # H (L1,L2)   #H(L1,L2)=((Linf-L1)/Linf-L2)^(M/2K)
  df.CA$H <- NA
  for(x2 in 1:(length(df.CA$H)-1)){
    df.CA$H[x2] <- ((Linf - df.CA$lowerLength.CA[x2]) /
                     (Linf - df.CA$lowerLength.CA[x2+1])) ^
      (M / (2*K))
  }

  #Survivors    #N(L1)=(N(L2)*H(L1,L2)+C(L1,L2)) *H(L1,L2)
  df.CA$survivors <- NA

  # survivors last size class
  df.CA$survivors[length(df.CA$survivors)] <-
    df.CA$catchCor.CA[length(df.CA$survivors)] / (terminalF/(terminalF + M))
  # other survivors
  for(x3 in (length(df.CA$survivors)-1):1){
    df.CA$survivors[x3] <- (df.CA$survivors[x3+1] *
                             df.CA$H[x3] + df.CA$catchCor.CA[x3] ) *
      df.CA$H[x3]
  }

  # F/Z  #F(L1,L2)/Z(L1,L2)=C(L1,L2)/(N(L1)-N(L2))
  df.CA$F_Z <- NA
  for(x4 in 1:(length(df.CA$F_Z)-1)){
    df.CA$F_Z[x4] <- df.CA$catchCor.CA[x4] /
      (df.CA$survivors[x4] - df.CA$survivors[x4+1])
  }
  df.CA$F_Z[length(df.CA$F_Z)] <- terminalF / (terminalF + M)

  #F  # F = M * (F_Z / 1-F_Z)
  df.CA$F <- NA
  for(x5 in 1:(length(df.CA$F))){
    df.CA$F[x5] <- M  *  (df.CA$F_Z[x5] / (1 - df.CA$F_Z[x5]))
  }

  # Z
  df.CA$Z <- NA
  for(x6 in 1:(length(df.CA$Z))){
    df.CA$Z[x6] <- M  +  df.CA$F[x6]
  }

  #Annual mean Nr
  df.CA$annualMeanNr <- NA
  for(x7 in 1:(length(df.CA$annualMeanNr-1))){
    df.CA$annualMeanNr[x7] <- (df.CA$survivors[x7] -
                                df.CA$survivors[x7+1]) / df.CA$Z[x7]
  }

  #Mean body weight
  df.CA$meanBodyWeight <- a * df.CA$midLengths.CA ^ b

  #Mean biomass
  df.CA$meanBiomass <- df.CA$annualMeanNr * df.CA$meanBodyWeight
  df.CA$meanBiomassTon <- df.CA$meanBiomass/1000

  #Yield
  df.CA$yield <- df.CA$catchCor.CA * df.CA$meanBodyWeight
  df.CA$yieldTon <- df.CA$yield/1000

  #FOR PLOT
  #Survivors rearranged
  df.CA$survivors_rea <- NA
  for(x8 in 1:(length(df.CA$survivors_rea)-1)){
    df.CA$survivors_rea[x8] <- df.CA$survivors[x8+1]
  }
  df.CA$survivors_rea[length(df.CA$survivors_rea)] <- 0

  #Calculate natural losses
  df.CA$natLoss <- NA
  for(x9 in 1:length(df.CA$natLoss)){
    df.CA$natLoss[x9] <- df.CA$survivors[x9] - df.CA$survivors_rea[x9] -
      df.CA$catchCor.CA[x9]
  }

  #put together in dataframe
  df.CAnew <- data.frame(survivors = df.CA$survivors_rea,
                         nat.losses = df.CA$natLoss,
                         catch = df.CA$catchCor.CA)

  #transpose matrix for barplot function
  df.CAnew <- t(as.matrix(df.CAnew))
  colnames(df.CAnew) <- df.CA$midLengths.CA

  #save x axis positions
  par(new = F)
  mids <- barplot(df.CAnew, xlab="",
                  ylim=c(0,(max(df.CA$survivors)+
                              max(df.CA$survivors)/14)))

  #create CA plot
  par(mar = c(5, 4, 4, 4) + 0.3)
  barplot(df.CAnew,col=c('darkgreen','purple','yellow'),
          ylim=c(0,(max(df.CA$survivors) + max(df.CA$survivors)/14)),
          axisnames = F,axis.lty = F,
          xlab = "Midlength [cm]", ylab = "Population" )                   #,names.arg = as.character(data$Midlength))#legend = rownames(data.new),
  legend('topright',legend = c(rownames(df.CAnew),"Fishing mortality"),
         fill = c('darkgreen','purple','yellow',NA),col='red',
         cex=0.9, xpd = T, text.width = 10,
         x.intersp = c(0.5,0.5,0.5,0.6),y.intersp = 1.4,
         lty = c(NA,NA,NA,1),lwd=3,merge=T,border = c(T,T,T,F),
         seg.len = 0.8)
  axis(1,at=mids[seq(1,length(mids),3)],
       labels=df.CA$midLengths.CA[seq(1,length(df.CA$midLengths.CA),3)],
       tick = F)
  par(new = TRUE)
  plot(df.CA$midLengths.CA, df.CA$F, type = "l", col='red',lwd=4,
       axes = FALSE, bty = "n", xlab = "", ylab = "")
  axis(side=4, at = pretty(range(df.CA$F)))
  mtext("Fishing mortatlity", side=4, line=3)
  plot1 <- recordPlot()

  #save all in list
  results.CA <- list()
  results.CA[[1]] <- df.CA
  results.CA[[2]] <- df.CAnew
  results.CA[[3]] <- plot1
  names(results.CA) <- c("Dataframe","Plotting_dataframe","Plot")

  return(results.CA)
}





