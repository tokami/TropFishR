#' @title Cumulative catch curve (The Jones and van Zalinge method)
#'
#' @description  This is a function to calculate the total mortality (Z) from length
#'    composition data via the length converted cumulative catch curve or from age
#'    composition data via the cumulative catch curve.
#'
#' @param classes Midpoints of the size class intervals as vector
#' @param cumCatch Cumulative catch per sampling time as matrix or the total cumulative catch as vector.
#' @param datatype Type of data which is used for analysis, either 'length' or 'age', for length frequency or age composition data, respectively
#' @param Linf Infinite length for investigated species in cm [cm].
#' @param K Growth coefficient for investigated species per year [1/year].
#' @param t0 Theoretical time zero, at which individuals of this species hatch (default = 0).
#'
#' @examples
#' \donttest{
#'
#' # Cumulative catch curve based on length frequency data
#' # load data
#'  data(ex.LengthCC)
#'
#' # calculate cumulative catch
#'  cumulativeCatch <- rev(cumsum(rev(ex.LengthCC[,2])))
#'
#' # run model
#'  output <- CatchCurveCum(classes = ex.LengthCC[,1],
#'    cumCatch = cumulativeCatch, datatype = "length", Linf = 23.1, K = 0.59)
#'
#' # investigate results
#'  output
#'
#' # based on age composition data
#' # calculate cumulative catch
#' cumulativeCatch <- rev(cumsum(rev(c(488,612,601,237,62.3,15.7,4.7,1.4))))
#'
#' # run model
#' output <- CatchCurveCum(classes = c(0,1,2,3,4,5,6,7), cumulativeCatch,
#'    datatype = 'age')
#'
#' # investigate results
#' output
#'
#' }
#'
#' @details Not good for methods where the largest length groups are not represented fully (highly selective gear). The length converted linearised cumulative catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @export

CatchCurveCum <- function(classes, cumCatch, datatype,
                          Linf = NULL, K = NULL, t0 = 0){

  # Error message if catch and age do not have same length
  if(class(cumCatch) == 'matrix' | class(cumCatch) == 'data.frame'){
    if(length(classes) != length(cumCatch[,1])) stop("Ages and catch do not have the same length!")
  }else if(class(cumCatch) == 'numeric'){
    if(length(classes) != length(cumCatch)) stop("Ages and catch do not have the same length!")
  }

  df.CCC <- cbind(classes,cumCatch)
  df.CCC <- as.data.frame(df.CCC)
  df.CCC$classes <- as.character(df.CCC$classes)

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(df.CCC$classes, split="\\+"))
  df.CCC$classes.num <- as.numeric(classes.num[,1])

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #    Linearised catch curve with variable time intervals   #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(class(cumCatch) == 'numeric'){

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #      Length converted catch curve     #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(datatype == 'length' & (is.null(Linf) | is.null(K))) stop("You need to assign values to Linf and K for the Catch curve based on length frequency data!")
    if(datatype == 'length'){

      #calculate size class interval
      df.CCC$midLengths <- df.CCC$classes.num
      interval.CCC <- df.CCC$midLengths[2] - df.CCC$midLengths[1]

      # L of lower length classes
      df.CCC$lowerLengths <- df.CCC$midLengths - (interval.CCC / 2)

      #ln C(L1,Linf)
      df.CCC$ln_C <- log(df.CCC$cumCatch)

      #ln (Linf - L)
      df.CCC$ln_Linf_L <- log(Linf - df.CCC$lowerLengths)

      #identify plot
      plot(x = df.CCC$ln_Linf_L,y = df.CCC$ln_C,
           xlab = "ln (Linf - L)", ylab = "ln C(L, Linf)")
      print(cat("Please choose the minimum and maximum point in the graph to include for the regression line!\nThen press 'Finish'!"))
      cutter <- identify(x = df.CCC$ln_Linf_L, y = df.CCC$ln_C,
                         labels = rownames(df.CCC),n=2)

      #calculations + model
      df.CCC.cut <- df.CCC[cutter[1]:cutter[2],]
      lm1 <- lm(ln_C ~ ln_Linf_L, data = df.CCC.cut)
      sum_lm1 <- summary(lm1)
      r_lm1 <- sum_lm1$r.squared
      intercept_lm1 <- sum_lm1$coefficients[1]
      slope_lm1 <- sum_lm1$coefficients[2]
      se_slope_lm1 <- sum_lm1$coefficients[4]

      #fit of regression line
      lm1.fit <- sum_lm1$r.squared

      ZK_lm1 <- abs(slope_lm1)
      SE_ZK_lm1 <- abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])

      Z_lm1 <- ZK_lm1 * K
      SE_Z_lm1 <- SE_ZK_lm1 * K

      #final plot
      plot(x = df.CCC$ln_Linf_L, y = df.CCC$ln_C,
           xlab = "ln (Linf - L)", ylab = "ln C(L, Linf)",cex = 1.5, bty = 'n')
      par(new=T)
      points(x = df.CCC.cut$ln_Linf_L, y = df.CCC.cut$ln_C,
             pch = 19, col = 'blue', cex = 1.5)
      segments(x0=df.CCC$ln_Linf_L[cutter[1]],
               y0=df.CCC$ln_C[cutter[1]],
               x1=df.CCC$ln_Linf_L[cutter[2]],
               y1=df.CCC$ln_C[cutter[2]],
               col="blue",lwd = 1.7)
      mtext(side = 3, text = paste("Z =",round(Z_lm1,2),"+/-",
                                   round(SE_Z_lm1,2)), col = 'blue')
      plot1 = recordPlot()

      #save all in list
      results.CCC <- list()
      results.CCC[[1]] <- df.CCC
      results.CCC[[2]] <- paste("Z =",round(Z_lm1,2),"+/-", round(SE_Z_lm1,2))
      results.CCC[[3]] <- plot1
      names(results.CCC) <- c("Dataframe","Total_mortality","Plot")

      return(results.CCC)
    }

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #    Aged based Catch curve    #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(datatype == 'age'){

      #NECESSARY???
#       # delta t
#       df.CCC$dt <- NA
#       for(x1 in 1:(length(df.CCC$dt)-1)){
#         df.CCC$dt[x1] <- df.CCC$classes.num[x1+1] - df.CCC$classes.num[x1]
#       }
#
#       # (t + dt) / 2   ==   x
#       df.CCC$tplusdt_2 <- (df.CCC$classes.num + df.CCC$dt ) / 2

      df.CCC$tplusdt_2 <- df.CCC$classes.num

      #ln C(L1,Linf)
      df.CCC$ln_C <- log(df.CCC$cumCatch)

      #identify plot
      plot(x = df.CCC$tplusdt_2,y = df.CCC$ln_C,
           xlab = "Age [yrs]", ylab = "ln C(t, inf)")
      print(cat("Please choose the minimum and maximum point in the graph to include for the regression line!\nThen press 'Finish'!"))
      cutter <- identify(x = df.CCC$tplusdt_2, y = df.CCC$ln_C,
                         labels = rownames(df.CCC),n=2)

      #calculations + model
      df.CCC.cut <- df.CCC[cutter[1]:cutter[2],]
      lm1 <- lm(ln_C ~ tplusdt_2, data = df.CCC.cut)
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
      plot(x = df.CCC$tplusdt_2, y = df.CCC$ln_C,
           xlab = "Age [yrs]", ylab = "ln C(t, inf)",cex = 1.5, bty = 'n')
      par(new=T)
      points(x = df.CCC.cut$tplusdt_2, y = df.CCC.cut$ln_C,
             pch = 19, col = 'blue', cex = 1.5)
      segments(x0=df.CCC$tplusdt_2[cutter[1]],
               y0=df.CCC$ln_C[cutter[1]],
               x1=df.CCC$tplusdt_2[cutter[2]],
               y1=df.CCC$ln_C[cutter[2]],
               col="blue",lwd = 1.7)
      mtext(side = 3, text = paste("Z =",round(Z_lm1,2),"+/-",
                                   round(SE_Z_lm1,2)), col = 'blue')
      plot1 = recordPlot()

      #save all in list
      results.CCC <- list()
      results.CCC[[1]] <- df.CCC
      results.CCC[[2]] <- paste("Z =",round(Z_lm1,2),"+/-", round(SE_Z_lm1,2))
      results.CCC[[3]] <- plot1
      names(results.CCC) <- c("Dataframe","Total_mortality","Plot")

      return(results.CCC)
    }
  }
}
