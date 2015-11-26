#' @title Catch curve
#
#' @description  This is a function to calculate the total mortality (Z) from length composition data via the length converted catch curve or from age at length data with catch curve.
#'
#' @param classes Midpoints of the length class as vector (length frequency data) or ages as vector (age composition data).
#' @param catch The catch as vector, or a matrix with catches of subsequent years if the catch curve with constat time intervals should be applied.
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficent for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch (default: 0).
#'
#' @examples
#' \donttest{
#' # Variable paramter system
#' # based on length frequency data
#' # load data
#' data("ex.LengthCC")
#'
#' # run model
#' output <- with(ex.LengthCC,CatchCurve(classes = ex.LengthCC[,1], catch = ex.LengthCC[,2],
#'   Linf = 23.1, K = 0.59, datatype = 'length'))
#'
#' # investigate results
#' output
#'
#' # based on age composition data
#' # load data
#' data("ex.CatchCurve")
#'
#' # run model
#' output <- with(ex.CatchCurve,CatchCurve(age, ex.CatchCurve[,2],datatype = 'age'))
#'
#' output
#'
#' # Constant parameter system based on age composition data
#' output <- with(ex.CatchCurve,CatchCurve(age, ex.CatchCurve[,2:8],datatype = 'age'))
#'
#' output
#'  }
#'
#' @details For variable parameter system vectors are reuqired for constant parameter systems matrices or data.frames have to be inserted. or vectors The length converted linearised catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references
#' Sparre ??? data from ICES 1981a
#'
#' @export

CatchCurve <- function(classes, catch, datatype, Linf = NULL, K = NULL, t0 = 0){

  # Error message if catch and age do not have same length
  if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
    if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
  }else if(class(catch) == 'numeric'){
    if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
  }

  df.CC <- cbind(classes,catch)
  df.CC <- as.data.frame(df.CC)
  df.CC$classes <- as.character(df.CC$classes)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(df.CC$classes, split="\\+"))
  df.CC$classes.num <- as.numeric(classes.num[,1])


  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #   Linearised catch curve with constant time intervals    #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(class(catch) == 'matrix' | class(catch) == 'data.frame'){

    if(datatype == 'length') stop("The catch curve with constant time interval is not applicable to length frequency data.")

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #   Aged based Catch curve   #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(datatype == 'age'){

      #find cohort to analyse
      real.cohort <- diag(as.matrix(catch))
      df.CC$catch.cohort <- c(real.cohort,
                              rep(NA,length(df.CC$classes.num) - length(real.cohort)))

      # ln( Catch )     ==    y
      df.CC$lnC <- log(df.CC$catch.cohort)

      #identify plot
      plot(x = df.CC$classes.num,y = df.CC$lnC,
           xlab = "Age [yrs]", ylab = "ln(C)")
      print("Please choose the minimum and maximum point in the graph to include for the regression line!")
      cutter <- identify(x = df.CC$classes.num, y = df.CC$lnC,
                         labels = rownames(df.CC), n=2)

      #calculations + model
      df.CC.cut <- df.CC[cutter[1]:cutter[2],]
      lm1 <- lm(lnC ~ classes.num, data = df.CC.cut)
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
      plot(x = df.CC$classes.num, y = df.CC$lnC,
           xlab = "Age [yrs]", ylab = "ln(C)",
           cex = 1.5, bty = 'n')
      par(new=T)
      points(x = df.CC.cut$classes.num, y = df.CC.cut$lnC,
             pch = 19, col = 'blue', cex = 1.5)
      segments(x0=df.CC$classes.num[cutter[1]],
               y0=df.CC$lnC[cutter[1]],
               x1=df.CC$classes.num[cutter[2]],
               y1=df.CC$lnC[cutter[2]],
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
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #    Linearised catch curve with variable time intervals   #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(class(catch) == 'numeric'){

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #      Length converted catch curve     #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(datatype == 'length' & (is.null(Linf) | is.null(K))) stop("You need to assign values to Linf and K for the Catch curve based on length frequency data!")
    if(datatype == 'length'){

      #calculate size class interval
      df.CC$midLengths.CC <- df.CC$classes.num
      interval.CC <- df.CC$midLengths.CC[2] - df.CC$midLengths.CC[1]

      # t of lower length classes
      df.CC$lowerLengths.CC <- df.CC$midLengths.CC - (interval.CC / 2)
      df.CC$t_L1 <- (t0 - (1/K)) * log(1 - (df.CC$lowerLengths.CC / Linf))

      # delta t
      df.CC$dt <- NA
      for(x1 in 1:(length(df.CC$dt)-1)){
        df.CC$dt[x1] <- df.CC$t_L1[x1+1] - df.CC$t_L1[x1]
      }

      # t of midlengths
      df.CC$t_midL <- (t0 - (1/K)) * log(1 - (df.CC$midLengths.CC / Linf))
      # ln( Catch / delta t)
      df.CC$lnC_dt <- log(df.CC$catch / df.CC$dt)


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
      SE_Z_lm1 <- abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])


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

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #    Aged based Catch curve    #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(datatype == 'age'){

      # delta t
      df.CC$dt <- NA
      for(x1 in 1:(length(df.CC$dt)-1)){
        df.CC$dt[x1] <- df.CC$classes.num[x1+1] - df.CC$classes.num[x1]
      }

      # (t + dt) / 2   ==   x
      df.CC$tplusdt_2 <- df.CC$classes.num + (df.CC$dt / 2)

      # ln( Catch / delta t)     ==    y
      df.CC$lnC_dt <- log(df.CC$catch / df.CC$dt)


      #identify plot
      plot(x = df.CC$tplusdt_2,y = df.CC$lnC_dt,
           xlab = "Age [yrs]", ylab = "ln(C/dt)")
      print("Please choose the minimum and maximum point in the graph to include for the regression line!")
      cutter <- identify(x = df.CC$tplusdt_2, y = df.CC$lnC_dt,
                         labels = rownames(df.CC), n=2)

      #calculations + model
      df.CC.cut <- df.CC[cutter[1]:cutter[2],]
      lm1 <- lm(lnC_dt ~ tplusdt_2, data = df.CC.cut)
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
      plot(x = df.CC$tplusdt_2, y = df.CC$lnC_dt,
           xlab = "Age [yrs]", ylab = "ln(C / dt)",
           cex = 1.5, bty = 'n')
      par(new=T)
      points(x = df.CC.cut$tplusdt_2, y = df.CC.cut$lnC_dt,
             pch = 19, col = 'blue', cex = 1.5)
      segments(x0=df.CC$tplusdt_2[cutter[1]],
               y0=df.CC$lnC_dt[cutter[1]],
               x1=df.CC$tplusdt_2[cutter[2]],
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
  }
}
