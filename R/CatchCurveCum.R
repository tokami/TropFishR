#' @title Cumulative catch curve (The Jones and van Zalinge method)
#'
#' @description  This is a function to calculate the total mortality (Z) from length
#'    composition data via the length converted cumulative catch curve or from age
#'    composition data via the cumulative catch curve.
#'
#' @param param A list consisting of following parameters:
#'   \code{$age} or \code{$midLengths} midpoints of the length class as vector (length frequency
#'   data) or ages as vector (age composition data),
#'   \code{$Linf} Infinite length for investigated species in cm [cm],
#'   \code{$K} Growth coefficent for investigated species per year [1/year],
#'   \code{t0} Theoretical time zero, at which individuals of this species hatch,
#'   \code{catch} Catch as vector, or a matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#' @param catch_column A number indicating the column of the catch matrix which should be
#'   used for the anaalysis.
#'
#' @examples
#' \donttest{
#'
#' # Cumulative catch curve based on length frequency data
#' # load data
#' data(goatfish)
#'
#' # calculate cumulative catch and add to data list
#' goatfish$cumCatch <- rev(cumsum(rev(goatfish$catch)))
#'
#' # run model
#' CatchCurveCum(goatfish)
#'
#' # based on age composition data
#' data(synCAA2)
#'
#' # calculate cumulative catch
#' synCAA2$cumCatch <- rev(cumsum(rev(synCAA2$catch)))
#'
#' # run model
#' CatchCurveCum(synCAA2)
#'
#' }
#'
#' @details Not good for methods where the largest length groups are not represented fully (highly selective gear). The length converted linearised cumulative catch curve is used to calculate the total mortality (Z). This function includes a so called locator function, which asks you to choose points from a graph manually. Based on these points the regression line is calculated.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

CatchCurveCum <- function(param){

  res <- param
  cumCatch <- res$cumCatch
  if(is.null(cumCatch)) stop("Please provide a cumulative catch!")
  if("midLengths" %in% names(res)) classes <- as.character(res$midLengths)
  if("age" %in% names(res)) classes <- as.character(res$age)


  # Error message if catch and age do not have same length
  if(class(cumCatch) == 'matrix' | class(cumCatch) == 'data.frame'){
    if(length(classes) != length(cumCatch[,1])) stop("Ages and catch do not have the same length!")
  }else if(class(cumCatch) == 'numeric'){
    if(length(classes) != length(cumCatch)) stop("Ages and catch do not have the same length!")
  }

  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #    Linearised catch curve with variable time intervals   #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if(class(cumCatch) == 'numeric'){

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #      Length converted catch curve     #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if("midLengths" %in% names(res) == TRUE){

      Linf <- res$Linf
      K <- res$K
      t0 <- ifelse("t0" %in% names(res),res$t0,0)

      if(is.null(Linf) | is.null(K)) stop("You need to assign values to Linf and K for the Catch curve based on length frequency data!")


      #calculate size class interval
      midLengths <- classes.num
      interval <- midLengths[2] - midLengths[1]

      # L of lower length classes
      lowerLengths <- midLengths - (interval / 2)

      #ln C(L1,Linf)
      ln_C <- log(cumCatch)

      #ln (Linf - L)
      ln_Linf_L <- log(Linf - lowerLengths)

      #for plot
      minlnC <- ifelse(min(lnC,na.rm=TRUE) < 0, min(lnC,na.rm=TRUE),0)
      maxlnC <- max(lnC,na.rm=TRUE) + 1

      #identify plot
      plot(x = ln_Linf_L,y = ln_C, ylim = c(minlnC,maxlnC),
           xlab = "ln (Linf - L)", ylab = "ln C(L, Linf)")
      print(cat("Please choose the minimum and maximum point in the graph to include for the regression line!\nThen press 'Finish'!"))
      cutter <- identify(x = ln_Linf_L, y = ln_C,
                         labels = order(ln_Linf_L),n=2)

      #calculations + model
      df.CCC <- as.data.frame(cbind(ln_Linf_L,ln_C))
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
      plot(x = ln_Linf_L, y = ln_C, ylim = c(minlnC,maxlnC),
           xlab = "ln (Linf - L)", ylab = "ln C(L, Linf)",cex = 1.5)
      par(new=T)
      points(x = df.CCC.cut$ln_Linf_L, y = df.CCC.cut$ln_C,
             pch = 19, col = 'blue', cex = 1.5)
      segments(x0=ln_Linf_L[cutter[1]],
               y0=ln_C[cutter[1]],
               x1=ln_Linf_L[cutter[2]],
               y1=ln_C[cutter[2]],
               col="blue",lwd = 1.7)
      mtext(side = 3, text = paste("Z =",round(Z_lm1,2),"+/-",
                                   round(SE_Z_lm1,2)), col = 'blue')

      #save all in list
      ret <- c(res,list(
        classes.num = classes.num,
        ln_Linf_L = ln_Linf_L,
        ln_C = ln_C,
        Z =  Z_lm1,
        se = SE_Z_lm1
      ))
      return(ret)
    }

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #    Aged based Catch curve    #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if("age" %in% names(res) == TRUE){

      #NECESSARY???
#       # delta t
#       dt <- NA
#       for(x1 in 1:(length(dt)-1)){
#         dt[x1] <- classes.num[x1+1] - classes.num[x1]
#       }
#
#       # (t + dt) / 2   ==   x
#       tplusdt_2 <- (classes.num + dt ) / 2

      tplusdt_2 <- classes.num

      #ln C(L1,Linf)
      ln_C <- log(cumCatch)

      #for plot
      minlnC <- ifelse(min(lnC,na.rm=TRUE) < 0, min(lnC,na.rm=TRUE),0)
      maxlnC <- max(lnC,na.rm=TRUE) + 1

      #identify plot
      plot(x = tplusdt_2,y = ln_C, ylim = c(minlnC,maxlnC),
           xlab = "Age [yrs]", ylab = "ln C(t, inf)")
      print(cat("Please choose the minimum and maximum point in the graph to include for the regression line!\nThen press 'Finish'!"))
      cutter <- identify(x = tplusdt_2, y = ln_C,
                         labels = order(tplusdt_2),n=2)

      #calculations + model
      df.CCC <- as.data.frame(cbind(tplusdt_2,ln_C))
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
      plot(x = tplusdt_2, y = ln_C, ylim = c(minlnC,maxlnC),
           xlab = "Age [yrs]", ylab = "ln C(t, inf)",cex = 1.5)
      par(new=T)
      points(x = df.CCC.cut$tplusdt_2, y = df.CCC.cut$ln_C,
             pch = 19, col = 'blue', cex = 1.5)
      segments(x0=tplusdt_2[cutter[1]],
               y0=ln_C[cutter[1]],
               x1=tplusdt_2[cutter[2]],
               y1=ln_C[cutter[2]],
               col="blue",lwd = 1.7)
      mtext(side = 3, text = paste("Z =",round(Z_lm1,2),"+/-",
                                   round(SE_Z_lm1,2)), col = 'blue')

      #save all in list
      ret <- c(res,list(
        tplusdt_2 = tplusdt_2,
        ln_C = ln_C,
        Z =  Z_lm1,
        se = SE_Z_lm1
      ))
      return(ret)
    }
  }
}
