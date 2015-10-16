#' @title Beverton and Holt's Z-Equations
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param classes Midpoints of the length class as vector (length frequency data) or ages as vector (age composition data).
#' @param catch Catch as vector, or a matrix with catches of subsequent years if the catch curve with constat time intervals should be applied.
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param PowellWetherall A logical parameter indicating if the Powell Wetherall method should be applied to length frequency data.
#'
#' @examples
#' data("ex.BevertonHoltsZ")
#' with(ex.BevertonHoltsZ, BevertonHoltsZ(midLength,
#' catch1960, datatype="length", Linf = 100, K = 0.3))
#' \donttest{
#' data("ex.PowellWetherall")
#' with(ex.PowellWetherall, BevertonHoltsZ(midLength,
#' catch, datatype="length", PowellWetherall = T))
#' data("ex.PowellWetherall2")
#' with(ex.PowellWetherall2, BevertonHoltsZ(midLength,
#' catch, datatype="length", PowellWetherall = T))
#' }
#'
#' @details For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references
#' Sparre ???
#'
#' @export

BevertonHoltsZ <- function(classes, catch, datatype, PowellWetherall = FALSE,
                           Linf = NULL, K = NULL){

  # Error message if catch and age do not have same length
  if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
    if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
  }else if(class(catch) == 'numeric'){
    if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
  }

  df.BH <- cbind(classes,catch)
  df.BH <- as.data.frame(df.BH)
  df.BH$classes <- as.character(df.BH$classes)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(df.BH$classes, split="\\+"))
  df.BH$classes.num <- as.numeric(classes.num[,1])

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #   Length based equation    #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(datatype == 'length'){
    if(PowellWetherall == FALSE){
      # calculate  C * (L1 + L2) / 2
      df.BH$c_midlength <- df.BH$catch * df.BH$classes.num

      # calculate L mean
      Lmean <- sum(df.BH$c_midlength) / sum(df.BH$catch)

      # calculate L prime
      Lprime <- df.BH$classes.num[1] -
        ((df.BH$classes.num[2] - df.BH$classes.num[1]) / 2 )

      Z.BH = K * (Linf - Lmean) / (Lmean - Lprime)

      #save all in list
      results.BH <- list()
      results.BH[[1]] <- df.BH
      results.BH[[2]] <- paste("Z =",round(Z.BH,2))
      names(results.BH) <- c("Dataframe","Total_mortality")

      return(results.BH)
    }


    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #    Powell-Wetherall method    #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(PowellWetherall == TRUE){

      # calculate cumulative catch
      df.BH$cumCatch <- rev(cumsum(rev(df.BH$catch)))

      # calculate  C * (L1 + L2) / 2
      df.BH$c_midlength <- df.BH$catch * df.BH$classes.num

      # calculate L prime   ===   x
      interval.BH <-  (df.BH$classes.num[2] - df.BH$classes.num[1]) / 2
      df.BH$Lprime <- df.BH$classes.num - interval.BH

      # calculate L mean
      for(i in 1:length(df.BH$c_midlength)){
        df.BH$sum_midL_c[i] <- sum(df.BH$c_midlength[i:length(df.BH$c_midlength)])
        df.BH$Lmean[i] <- sum(df.BH$c_midlength[i:length(df.BH$c_midlength)]) /
          sum(df.BH$catch[i:length(df.BH$catch)])
      }

      # calculate Lmean - Lprime
      df.BH$Lmean_Lprime <- df.BH$Lmean - df.BH$Lprime



      #identify plot
      plot(x = df.BH$Lprime,y = df.BH$Lmean_Lprime,
           xlab = "Lprime", ylab = "Lmean - Lprime")
      print("Please choose the minimum and maximum point in the graph to include for the regression line!")
      cutter <- identify(x = df.BH$Lprime,y = df.BH$Lmean_Lprime,
                         labels = rownames(df.BH), n=2)

      #calculations + model
      df.BH.cut <- df.BH[cutter[1]:cutter[2],]
      lm1 <- lm(Lmean_Lprime ~ Lprime, data = df.BH.cut)
      sum_lm1 <- summary(lm1)
      r_lm1 <- sum_lm1$r.squared
      intercept_lm1 <- sum_lm1$coefficients[1]
      slope_lm1 <- sum_lm1$coefficients[2]
      se_slope_lm1 <- sum_lm1$coefficients[4]
      se_intercept_lm1 <- sum_lm1$coefficients[3]

      #fit of regression line
      lm1.fit <- sum_lm1$r.squared

      SE_slope <- abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])
      SE_intercept <- abs(se_intercept_lm1) * qt(0.975,sum_lm1$df[1])

      Linf.BH <- - intercept_lm1 / slope_lm1
      ZK.BH <- - (1+slope_lm1)/slope_lm1

      #final plot
      plot(x = df.BH$Lprime,y = df.BH$Lmean_Lprime,
           xlab = "Lprime", ylab = "Lmean- Lprime",
           cex = 1.5, bty = 'n')
      par(new=T)
      points(x = df.BH.cut$Lprime,y = df.BH.cut$Lmean_Lprime,
             pch = 19, col = 'blue', cex = 1.5)
      abline(a=intercept_lm1,b=slope_lm1,col="blue",lwd = 1.7)
      #mtext(side = 3, text = paste("Z =",round(Z_lm1,2),"±",
      #                             round(SE_Z_lm1,2)), col = 'blue')
      plot1 = recordPlot()

      #save all in list
      results.BH <- list()
      results.BH[[1]] <- df.BH
      results.BH[[2]] <- paste("Z/K =",round(ZK.BH,2))
      results.BH[[3]] <- paste("Linf =",round(Linf.BH,2))
      results.BH[[4]] <- plot1
      names(results.BH) <- c("Dataframe","Total_mortality_K","Linf","Plot")

      return(results.BH)
    }
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #     Aged based equation    #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(datatype == 'age'){}

}
