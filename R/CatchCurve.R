#' @title Catch curve
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
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
#'  data("ex.CatchCurve")
#'  output <- with(ex.CatchCurve,
#'  LengthConCatchCurve(age, ex.CatchCurve[,2:8])
#'  output
#'  }
#'
#' @details The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references
#' Sparre ??? data from ICES 1981a
#'
#' @export

data("ex.CatchCurve")

age = ex.CatchCurve$age
catch = ex.CatchCurve[,2:8]

CatchCurve <- function(age, catch){


  if(length(age) != length(catch[,1])) stop("Ages and catch do not have the same length!")


  df.CC <- cbind(age,catch)
  df.CC$age <- as.character(df.CC$age)
  age.num <- do.call(rbind,strsplit(df.CC$age, split='+'))
  age.num <- as.numeric(age.num[,1])

  delta.age <- NA
  for(i in 1:(length(age.num)-1)){
    delta.age[i] <- age.num[i+1] - age.num[i]
  }


  if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
    for(i in 1:length(catch)){

      ln_C <- log(catch[i] / delta.age)

      #identify plot
      plot(x = df.CC$age.num, y = ln_C,
           xlab = "Relative age [yrs]", ylab = "ln(C/dt)")
      print("Please choose the minimum and maximum point in the graph to include for the regression line!")
      cutter <- identify(x = df.CC$t_midL, y = df.CC$lnC_dt,
                         labels = rownames(df.CC), n=2)

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

    }
  }else

  #identify plot
  plot(x = df.CC$t_midL,y = df.CC$lnC_dt,
       xlab = "Relative age [yrs]", ylab = "ln(C/dt)")
  print("Please choose the minimum and maximum point in the graph to include for the regression line!")
  cutter <- identify(x = df.CC$t_midL, y = df.CC$lnC_dt,
                     labels = rownames(df.CC), n=2)

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
