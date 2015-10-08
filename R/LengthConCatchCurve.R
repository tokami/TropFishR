#' @title Length converted catch curve
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve.
#'
#' @param midLengths Midpoints of the length class as vector
#' @param catch Catch per sampling time as matrix or the total catch as vector.
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch.
#' @param catchCorFac optional: Correction factor for catch, in case provided catch does spatially or temporarily not reflect catch for fishing ground of a whole year.
#'
#' @examples
#' \donttest{
#'  data("ex.LengthCC")
#'  output <- LengthConCatchCurve(midLengths = ex.LengthCC[,1], catch = ex.LengthCC[,2], Linf = 23.1, K = 0.59)
#'  output
#'  }
#'
#' @details The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references
#' Sparre ???
#'
#' @export

LengthConCatchCurve <- function(midLengths, catch, Linf, K, t0 = 0,
                       catchCorFac = NA){

  #Transform matrix into vector if provided
  if(class(catch) == 'matrix'){
    catch.vec <- rowSums(catch)
  }else catch.vec = catch

  if(length(midLengths) != length(catch.vec)) stop("midLengths and catch do not have the same length!")

  df.CC <- data.frame(midLengths.CC = midLengths,
                       catch.CC = catch.vec)

  #calculate size class interval
  interval.CC <- df.CC$midLengths.CC[2] - df.CC$midLengths.CC[1]

  # t of lower length classes
  df.CC$lowerLengths.CC <- df.CC$midLengths.CC - (interval.CC / 2)

  ### NECESSARY???
  if(!is.na(catchCorFac)){
    df.CC$catchCor.CC <- df.CC$catch.CC * catchCorFac
    }else df.CC$catchCor.CC <- df.CC$catch.CC

  df.CC$t_L1 <- (t0 - (1/K)) * log(1 - (df.CC$lowerLengths.CC / Linf))

  # delta t
  df.CC$dt <- NA
  for(x1 in 1:(length(df.CC$dt)-1)){
    df.CC$dt[x1] <- df.CC$t_L1[x1+1] - df.CC$t_L1[x1]
  }

  # t of midlengths
  df.CC$t_midL <- (t0 - (1/K)) * log(1 - (df.CC$midLengths.CC / Linf))
  # ln( Catch / delta t)
  df.CC$lnC_dt <- log(df.CC$catchCor.CC / df.CC$dt)


  #identify plot
  plot(x = df.CC$t_midL,y = df.CC$lnC_dt,
       xlab = "Relative age [yrs]", ylab = "ln(C/dt)")
  print("Please choose the minimum and maximum point in the graph to include for the regression line! Then press 'Finish'!")
  cutter <- identify(x = df.CC$t_midL, y = df.CC$lnC_dt,
                     labels = rownames(df.CC))

  #calculations + model
  df.CC.cut <- df.CC[cutter[1]:cutter[2],]
  lm1 <- lm(lnC_dt ~ t_midL, data = df.CC.cut)
  sum_lm1 <- summary(lm1)
  r_lm1 <- sum_lm1$r.squared
  intercept_lm1 <- sum_lm1$coefficients[1]
  slope_lm1 <- sum_lm1$coefficients[2]
  se_slope_lm1 <- sum_lm1$coefficients[4]

  #fit of regression line
  lm1.fit <- sum_lm1$r.squared

  Z_lm1 <- abs(slope_lm1)
  SE_Z_lm1 <- abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])   ##COMPARE VALUES WITH: sb2 <- ((sd.Y/sd.X)^2 - slope^2 )/ (length(n)-2)


  #final plot
  plot(x = df.CC$t_midL, y = df.CC$lnC_dt,
       xlab = "Relative age [yrs]", ylab = "ln(C / dt)",
       cex = 1.5, bty = 'n')
  par(new=T)
  points(x = df.CC.cut$t_midL, y = df.CC.cut$lnC_dt,
         pch = 19, col = 'blue', cex = 1.5)
  segments(x0=df.CC$t_midL[cutter[1]],
           y0=df.CC$lnC_dt[cutter[1]],
           x1=df.CC$t_midL[cutter[2]],
           y1=df.CC$lnC_dt[cutter[2]],
           col="blue",lwd = 1.7)
  mtext(side = 3, text = paste("Z =",round(Z_lm1,2),"±",
                               round(SE_Z_lm1,2)), col = 'blue')
  plot1 = recordPlot()


  #save all in list
  results.CC <- list()
  results.CC[[1]] <- df.CC
  results.CC[[2]] <- paste("Z =",round(Z_lm1,2),"±", round(SE_Z_lm1,2))
  results.CC[[3]] <- plot1
  names(results.CC) <- c("Dataframe","Total_mortality","Plot")

  return(results.CC)
}
