#' @title Beverton and Holt's Z-Equations & Powell and Wetherall method
#
#' @description A method to calculate Z from a range of equations derived by Beverton and
#'    Holt (1956) with the option to perform the Powell and Wetherall method (Powell, 1979;
#'    Wetherall et al., 1987). Refers to Chapter 4.5 in Sparre & Venema (1998).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{age} or \strong{midLengths}: midpoints of the length class as vector (length-frequency
#'   data) or ages as vector (age composition data),
#'   \item \strong{Linf}: infinite length for investigated species in cm [cm],
#'   \item \strong{K}: growth coefficent for investigated species per year [1/year],
#'   \item \strong{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \strong{catch}: catch as vector, or a matrix with catches of subsequent years;
#' }
#' @param catch_column optional; in case catch is a matrix or data.frame, a number
#'    indicating which column of the matrix should be analysed (Default: \code{NA}).
#' @param PowellWetherall logical; if \code{TRUE} the Powell Wetherall method is
#'   applied (Default: \code{FALSE}).
#'
#' @keywords function mortality Z
#'
#' @examples
#' #_______________________________________________
#' # Beverton and Holt's Z method
#' # based on length data
#' # load data
#' data(synLFQ2)
#'
#' # run model
#' Z_BevertonHolt(synLFQ2, catch_column = 2)
#'
#' # based on age data
#' # load data
#' data(synCAA1)
#'
#' #run model
#' Z_BevertonHolt(synCAA1, catch_column = 3)
#'
#' \donttest{
#' #_______________________________________________
#' # Powell and Wetherall method
#' # load data
#' data(synLFQ3)
#'
#' # run model
#' Z_BevertonHolt(synLFQ3, PowellWetherall = T)
#' }
#'
#' @details  The first length group or age class within the list object \strong{midLengths} or
#'    \strong{age} will be used as the Lprim or tprime (length of recruitment to fishery).
#'    In the case of the Powell and Wetherall method, this function includes a so called
#'    'locator' function, which asks you to choose two points from a graph manually. The
#'    two points which you choose by clicking on the plot in the graphical device represent
#'    the start and end of the data points, which should be used for the analysis. Based
#'    on these points the regression line is calculated. The Powell and Wetherall method
#'    only works with length-frequency data.
#'
#' @return A list with the input parameters and follwing objects:
#' \itemize{
#'   \item \strong{tmean} or \strong{Lmean}: mean age or length of fish,
#'   \item \strong{tprime} or \strong{Lprime}: some age or length for which all fish of that length and
#'      longer are under full exploitation,
#'   \item \strong{Z}: total mortality;}
#' and/or following objects when applying the Powell and Wetherall method:
#' \itemize{
#'   \item \strong{Lmean_Lprime}: dependent variable for regression analysis,
#'   \item \strong{Linf_calc}: infinite length in [cm],
#'   \item \strong{Z_K}: total mortality divided by K;
#' }
#'
#' @references
#' Beverton R.J.H and S.J. Holt, 1956. A review of methods of estimating mortality rates
#' in exploited fish populations, with special reference to sources of bias in catch
#' sampling. \emph{Rapp.P.-v.Réun.CIEM}, 140:67-83
#'
#' Powell, D.G., 1979. Estimation of mortality and growth parameters from the length-
#' frequency of a catch. \emph{Rapp.P.-v.Réun.CIEM}, 175:167-169
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' Wetherall, J.A., J.J. Polovina and S. Ralston, 1987. Estimating growth and mortality
#' in steady-state fish stocks from length-frequency data. \emph{ICLARM Conf.Proc.},
#' (13):53-74
#'
#' @export

Z_BevertonHolt <- function(param, catch_column = NA, PowellWetherall = FALSE){

  res <- param
  catch <- res$catch

  if(class(catch) == "data.frame" | class(catch) == "matrix"){
    if(is.na(catch_column)) stop("Please provide a number indicating which column of the catch matrix should be analysed!")
    catch <- catch[,catch_column]
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #   Length based equation    #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if("midLengths" %in% names(res)){

    classes <- as.character(res$midLengths)
    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    Linf <- res$Linf
    K <- res$K

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
    }

    if(PowellWetherall == FALSE){
      # calculate  C * (L1 + L2) / 2
      c_midlength <- catch * classes.num

      # calculate L mean
      Lmean <- sum(c_midlength) / sum(catch)

      # calculate L prime
      Lprime <- classes.num[1] -
        ((classes.num[2] - classes.num[1]) / 2 )

      Z = K * (Linf - Lmean) / (Lmean - Lprime)

      #save all in list
      ret <- c(res,list(
        Lmean = Lmean,
        Lprime = Lprime,
        Z = Z
      ))
      return(ret)
    }


    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #    Powell-Wetherall method    #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(PowellWetherall == TRUE){

      # calculate cumulative catch
      cumCatch <- rev(cumsum(rev(catch)))

      # calculate  C * (L1 + L2) / 2
      c_midlength <- catch * classes.num

      # calculate L prime   ===   x
      interval <-  (classes.num[2] - classes.num[1]) / 2
      Lprime <- classes.num - interval

      # calculate L mean
      sum_midL_c <- rep(NA,length(classes.num))
      Lmean <- rep(NA,length(classes.num))
      for(i in 1:length(c_midlength)){
        sum_midL_c[i] <- sum(c_midlength[i:length(c_midlength)])
        Lmean[i] <- sum(c_midlength[i:length(c_midlength)]) /
          sum(catch[i:length(catch)])
      }

      # calculate Lmean - Lprime
      Lmean_Lprime <- Lmean - Lprime

      #identify plot
      plot(x = Lprime,y = Lmean_Lprime,
           xlab = "Lprime", ylab = "Lmean - Lprime")
      print("Please choose the minimum and maximum point in the graph to include for the regression line!")
      cutter <- identify(x = Lprime,y = Lmean_Lprime,
                         labels = order(Lprime), n=2)

      #calculations + model
      df.BH <- as.data.frame(cbind(classes.num,Lmean_Lprime,Lprime))
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
      plot(x = Lprime,y = Lmean_Lprime,
           xlab = "Lprime", ylab = "Lmean- Lprime",
           cex = 1.5, main = "Powell-Wetherall plot")
      par(new=T)
      points(x = df.BH.cut$Lprime,y = df.BH.cut$Lmean_Lprime,
             pch = 19, col = 'blue', cex = 1.5)
      abline(a=intercept_lm1,b=slope_lm1,col="blue",lwd = 1.7)

      #save all in list
      ret <- c(res,list(
        Lmean_Lprime = Lmean_Lprime,
        Lprime = Lprime,
        Linf_calc = Linf.BH,
        Z_K = ZK.BH
      ))
      return(ret)
    }
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #     Aged based equation    #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if("midAge" %in% names(res) | "age" %in% names(res)){

    if("midAge" %in% names(res)) classes <- as.character(res$midAge)
    if("age" %in% names(res)) classes <- as.character(res$age)
    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
    }

    sample.size <- sum(catch,na.rm=T)
    sum.age.number <- sum((catch * classes.num), na.rm=T)
    tmean <- sum.age.number/sample.size
    interval <- (classes.num[2] - classes.num[1]) / 2
    tprime <- classes.num[1] - interval

    Z.BH <- 1 / (tmean - tprime)

    #save all in list
    ret <- c(res,list(
      tmean = tmean,
      tprime = tprime,
      Z = Z.BH
    ))
    return(ret)
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #     Length at first capture data    #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
}
