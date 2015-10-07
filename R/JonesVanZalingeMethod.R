#' @title Length converted cumulative catch curve (The Jones and van Zalinge method)
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted cumulative catch curve.
#'
#' @param midLengths Midpoints of the size class intervals as vector
#'
#' @param cumCatch Cumulative catch per sampling time as matrix or the total cumulative catch as vector.
#'
#' @param Linf Infinite length for investigated species in cm [cm].
#'
#' @param K Growth coefficient for investigated species per year [1/year].
#'
#' @param t0 Theoretical time zero, at which individuals of this species hatch.
#'
#' @param catchCorFac optional: Correction factor for catch, in case provided catch does spatially or temporarily not reflect catch for fishing ground of a whole year.
#'
#' @examples
#' \donttest{
#'  data("ex.LengthCC")
#'  #calculating the cumulated catch of fish of length L and above:
#'  cumulativeCatch <- rev(cumsum(rev(ex.LengthCC[,2])))
#'  output <- JonesVanZalingeMethod(midLengths = ex.LengthCC[,1],
#'    cumCatch = cumulativeCatch, Linf = 23.1, K = 0.59)
#'  output
#'  }
#'
#' @details The length converted linearised cumulative catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @export


JonesVanZalingeMethod <- function(midLengths, cumCatch, Linf, K, t0 = 0,
                                catchCorFac = NA){

  #Transform matrix into vector if provided
  if(class(cumCatch) == 'matrix'){
    cumCatch.vec <- rowSums(cumCatch)
  }else cumCatch.vec = cumCatch

  if(length(midLengths) != length(cumCatch.vec)) stop("midLengths and cumCatch do not have the same length!")

  df.JZ <- data.frame(midLengths.JZ = midLengths,
                      cumCatch.JZ = cumCatch.vec)

  #calculate size class interval
  interval.JZ <- df.JZ$midLengths.JZ[2] - df.JZ$midLengths.JZ[1]

  # L of lower length classes
  df.JZ$lowerLengths.JZ <- df.JZ$midLengths.JZ - (interval.JZ / 2)

  #ln C(L1,Linf)
  df.JZ$ln_C <- log(df.JZ$cumCatch.JZ)

  #ln (Linf - L)
  df.JZ$ln_Linf_L <- log(Linf - df.JZ$lowerLengths.JZ)

  ### NECESSARY???
  if(!is.na(catchCorFac)){
    df.JZ$cumCatchCor.JZ <- df.JZ$cumCatch.JZ * catchCorFac
  }else df.JZ$cumCatchCor.JZ <- df.JZ$cumCatch.JZ

  #identify plot
  plot(x = df.JZ$ln_Linf_L,y = df.JZ$ln_C,
       xlab = "ln (Linf - L)", ylab = "ln C(L, Linf)")
  print("Please choose the minimum and maximum point in the graph to include for the regression line! Then press 'Finish'!")
  cutter <- identify(x = df.JZ$ln_Linf_L, y = df.JZ$ln_C,
                     labels = rownames(df.JZ))

  #calculations + model
  df.JZ.cut <- df.JZ[cutter[1]:cutter[2],]
  lm1 <- lm(ln_C ~ ln_Linf_L, data = df.JZ.cut)
  sum_lm1 <- summary(lm1)
  r_lm1 <- sum_lm1$r.squared
  intercept_lm1 <- sum_lm1$coefficients[1]
  slope_lm1 <- sum_lm1$coefficients[2]
  se_slope_lm1 <- sum_lm1$coefficients[4]

  #fit of regression line
  lm1.fit <- sum_lm1$r.squared

  Z_lm1 <- abs(slope_lm1)
  SE_Z_lm1 <- abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])


  #final plot
  plot(x = df.JZ$ln_Linf_L, y = df.JZ$ln_C,
       xlab = "ln (Linf - L)", ylab = "ln C(L, Linf)",cex = 1.5, bty = 'n')
  par(new=T)
  points(x = df.JZ.cut$ln_Linf_L, y = df.JZ.cut$ln_C,
         pch = 19, col = 'blue', cex = 1.5)
  segments(x0=df.JZ$ln_Linf_L[cutter[1]],
           y0=df.JZ$ln_C[cutter[1]],
           x1=df.JZ$ln_Linf_L[cutter[2]],
           y1=df.JZ$ln_C[cutter[2]],
           col="blue",lwd = 1.7)
  mtext(side = 3, text = paste("Z =",round(Z_lm1,2),"±",
                               round(SE_Z_lm1,2)), col = 'blue')
  plot1 = recordPlot()

  #save all in list
  results.JZ <- list()
  results.JZ[[1]] <- df.JZ
  results.JZ[[2]] <- paste("Z =",round(Z_lm1,2),"±", round(SE_Z_lm1,2))
  results.JZ[[3]] <- plot1
  names(results.JZ) <- c("Dataframe","Total_mortality","Plot")

  return(results.JZ)
}
