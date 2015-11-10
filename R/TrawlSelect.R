#' @title Trawl net selctivity
#
#' @description  This function estimates the selecitvity of trawl nets.
#'
#' @param classes Midpoints of the length class as vector (length frequency data) or ages as vector (age composition data).
#' @param catch Catch as vector, or a matrix with catches of subsequent years if the catch curve with constat time intervals should be applied.
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch (default: 0).
#'
#' @examples
#' data("ex.TrawlSelect")
#' output <- with(ex.TrawlSelect,TrawlSelect(classes = midLengths,
#'    numCodend, numCover, meshsizeCodend = 4))
#' output
#'
#' @details To calculate selection factor (SF), L25, L50 and L75 for trawl nets /fisheries.
#'
#' @references xx
#'
#'
#' @export

TrawlSelect <- function(classes, numCodend, numCover, meshsizeCodend){

  df.TS <- cbind(classes,numCodend,numCover)
  df.TS <- as.data.frame(df.TS)
  df.TS$classes <- as.character(df.TS$classes)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(df.TS$classes, split="\\+"))
  df.TS$classes.num <- as.numeric(classes.num[,1])

  # calculate total number in codend and cover
  df.TS$numTot <- df.TS$numCodend + df.TS$numCover

  # calculate fraction retained (SL obs)
  df.TS$SLobs <- df.TS$numCodend/df.TS$numTot

  # ln( 1 / SL - 1)
  df.TS$lnSL <- log(1/df.TS$SLobs - 1)

  #excluding point where no (0) or full (1) retention was obtained and all values beyond those points, even if they are between 0 and 1
  df.TS$lnSL[which(df.TS$classes.num >= df.TS$classes.num[which(df.TS$SLobs == 1)])] <- NA
  df.TS$lnSL[which(df.TS$lnSL == Inf | df.TS$lnSL == -Inf)] <- NA



  #model
  mod.TS1 <- lm(lnSL ~ classes.num,df.TS)
  sum.mod.TS1 <- summary(mod.TS1)
  S1.TS <- sum.mod.TS1$coefficients[1]
  S2.TS <- abs(sum.mod.TS1$coefficients[2])

  #plot
  plot(df.TS$classes.num,df.TS$lnSL, main = "Regression analysis")
  abline(a=sum.mod.TS1$coefficients[1], b=sum.mod.TS1$coefficients[2])
  plot1 = recordPlot()

  #L25,L50,L75
  L25.TS <- (S1.TS - log(3)) / S2.TS
  L50.TS <- S1.TS / S2.TS
  L75.TS <- (S1.TS + log(3)) / S2.TS

  #estimated SL (SL est)
  df.TS$SLest <- 1 / (1 + exp(S1.TS - S2.TS * df.TS$classes.num))

  #Selection factor
  SF.TS <- L50.TS/meshsizeCodend

  #new plot
  xL <- seq(min(df.TS$classes.num,na.rm=T),max(df.TS$classes.num,na.rm=T),0.1)
  plot(SLobs ~ classes.num, df.TS, main = "Gear selection ogive",
       xlab = "L [cm]", ylab = "fraction retained", pch = 16 , cex =1.4)
  lines(x = xL, y = 1/(1 + exp(S1.TS - S2.TS * xL)), col='blue', lwd = 1.5)
  segments(x0 = L25.TS, x1 = L25.TS, y0 = 0,
           y1 = 1/(1 + exp(S1.TS - S2.TS * L25.TS)),col= 'gray40',lty = 2, lwd=1.5)
  segments(x0 = 0, x1 = L25.TS, y0 = 1/(1 + exp(S1.TS - S2.TS * L25.TS)),
           y1 = 1/(1 + exp(S1.TS - S2.TS * L25.TS)),col= 'gray40',lty = 2, lwd=1.5)
  text(labels = "L25%", x = L25.TS, y = (1/(1 + exp(S1.TS - S2.TS * L25.TS)))/3,
       col = 'gray40')
  segments(x0 = L50.TS, x1 = L50.TS, y0 = 0,
           y1 = 1/(1 + exp(S1.TS - S2.TS * L50.TS)),col= 'gray40',lty = 2, lwd=1.5)
  segments(x0 = 0, x1 = L50.TS, y0 = 1/(1 + exp(S1.TS - S2.TS * L50.TS)),
           y1 = 1/(1 + exp(S1.TS - S2.TS * L50.TS)),col= 'gray40',lty = 2, lwd=1.5)
  text(labels = "L50%", x = L50.TS, y = (1/(1 + exp(S1.TS - S2.TS * L50.TS)))/3,
       col = 'gray40')
  segments(x0 = L75.TS, x1 = L75.TS, y0 = 0,
           y1 = 1/(1 + exp(S1.TS - S2.TS * L75.TS)),col= 'gray40',lty = 2, lwd=1.5)
  segments(x0 = 0, x1 = L75.TS, y0 = 1/(1 + exp(S1.TS - S2.TS * L75.TS)),
           y1 = 1/(1 + exp(S1.TS - S2.TS * L75.TS)),col= 'gray40',lty = 2, lwd=1.5)
  text(labels = "L75%", x = L75.TS, y = (1/(1 + exp(S1.TS - S2.TS * L75.TS)))/3,
       col = 'gray40')

  plot2 = recordPlot()

  #all results together
  df.results <- data.frame(
    "S1" = S1.TS,
    "S2" = S2.TS,
    "L25" = L25.TS,
    "L50" = L50.TS,
    "L75" = L75.TS,
    "SF" = SF.TS
  )

  #save all in list
  results.TS <- list()
  results.TS[[1]] <- df.TS
  results.TS[[2]] <- df.results
  results.TS[[3]] <- plot1
  results.TS[[4]] <- plot2

  names(results.TS) <- c("Dataframe","Selectivity_Results",
                         "Regression_plot","Gear_selection_ogive")

  return(results.TS)
}
