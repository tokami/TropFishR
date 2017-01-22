#' @title Catch curve
#'
#' @description  This function applies the (length-converted) linearised catch
#'    curve to age composition and length-frequency data,
#'    respectively. It allows to estimate the instantaneous total mortality rate (Z).
#'    Optionally, the gear selectivity can be estimated and the cumulative catch
#'    curve cna be applied.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \code{midLengths} or \code{age}: midpoints of the length classes (length-frequency
#'   data) or ages (age composition data),
#'   \item \code{Linf}: infinite length for investigated species in cm [cm],
#'   \item \code{K}: growth coefficent for investigated species per year [1/year],
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \code{catch}: catches, vector or matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#' }
#' @param catch_columns numerical; indicating the column of the catch matrix which should be
#'   used for the analysis.
#' @param cumulative logical; if TRUE the cumulative
#'   catch curve is applied (Jones and van Zalinge method)
#' @param calc_ogive logical; if TRUE the selection ogive is additionally
#'   calculated from the catch curve (only if \code{cumulative = FALSE})
#' @param reg_int instead of using the identity method a range can be determined,
#'    which is to be used for the regression analysis. If equal to NULL identity method
#'    is applied (default).
#'
#' @keywords function mortality Z catchCurve
#'
#' @examples
#' \donttest{
#' #_______________________________________________
#' # Variable paramter system (with catch vector)
#' # based on length frequency data
#' data(goatfish)
#' output <- catchCurve(goatfish)
#' summary(output$linear_mod)
#'
#' # based on age composition data
#' data(whiting)
#' catchCurve(whiting, catch_columns = 1)
#'
#' #_______________________________________________
#' # Constant parameter system based on age composition data (with catch matrix)
#' catchCurve(whiting)
#'
#' #_______________________________________________
#' # Cumulative Catch Curve
#' # based on length frequency data
#' data(goatfish)
#' catchCurve(goatfish, cumulative = TRUE)
#'
#' # based on age composition data
#' data(synCAA2)
#' catchCurve(synCAA2, cumulative = TRUE)
#'
#' #_______________________________________________
#' # Catch Curve with estimation of selection ogive
#' data(synLFQ3)
#' output <- catchCurve(synLFQ3, calc_ogive = TRUE)
#' summary(output$linear_mod_sel)
#'  }
#'
#' # the same with predefined selection for regression line:
#' output <- catchCurve(synLFQ3, calc_ogive = TRUE, reg_int = c(9,21))
#' plot(output, plot_selec = TRUE)
#'
#' @details This function includes the \link{identify} function, which asks you to
#'   choose two points from a graph manually. The two points which you choose by clicking
#'   on the plot in the graphical device represent the start and end of the data points,
#'   which should be used for the analysis. Based on these points the regression line
#'   is calculated.
#'   When the selection ogive
#'   is calculated by means of the catch curve the assumption is made, that Z is constant
#'   for all year classes or length groups, respectively. Accoring to Sparre and Venema
#'   (1998) this assumption might be true, because F is smaller for young fish
#'   (Selectivity) while M is higher for young fish (high natural mortality). The selectivity
#'   for not fully exploited old fish (e.g. due to gillnet fishery) can not be calculated yet
#'   by use of the catch curve.
#'   Based on the format of the list argument \code{catch} and whether the argument
#'   \code{catch_columns} is defined, the function automatically
#'   distinguishes between the catch curve with variable parameter system (if catch is a
#'   vector) and the one with constant parameter system (if catch is a matrix or a
#'   data.frame and \code{catch_columns = NA}). In the case of the variable parameter
#'   system the catches of one year are
#'   assumed to represent the catches during the entire life span of a so called
#'   pseudo-cohort.
#'   The cumulative catch curve does not allow for the estimation of the selectivity
#'   ogive.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{classes.num}, \strong{tplusdt_2}, \strong{t_midL}, or
#'      \strong{ln_Linf_L}: age, relative age or subsitute depending on input and method,
#'   \item \strong{lnC} or \strong{lnC_dt}: logarithm of (rearranged) catches,
#'   \item \strong{reg_int}: the interval used for the regression analysis,
#'   \item \strong{linear_mod}: linear model used for the regression analysis,
#'   \item \strong{Z}: instantaneous total mortality rate, confidenceInt
#'   \item \strong{se}: standard error of the total mortality;
#'   \item \strong{confidenceInt}: confidence interval of the total mortality;}
#' in case calc_ogive == TRUE, additionally:
#' \itemize{
#'   \item \strong{intercept}: intercept of regression analysis,
#'   \item \strong{linear_mod_sel}: linear model used for the selectivity analysis,
#'   \item \strong{Sobs}: observed selection ogive,
#'   \item \strong{ln_1_S_1}: dependent variable of regression analysis for
#'   selectivity parameters,
#'   \item \strong{Sest}: estimated selection ogive,
#'   \item \strong{t50}: age at first capture (age at which fish have a 50%
#'   probability to be caught),
#'   \item \strong{t75}: age at which fish have a 75% probability to be caught,
#'   \item \strong{L50}: length at first capture (length at which fish have a 50%
#'   probability to be caught),
#'   \item \strong{L75}: length at which fish have a 75% probability to be caught;
#' }
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics identify par plot
#' @importFrom stats lm na.omit
#' @importFrom utils flush.console
#'
#' @references
#' Baranov, F.I., 1926. On the question of the dynamics of the fishing industry.
#' \emph{Nauchn. Byull. Rybn. Khoz}, 8 (1925), 7-11
#'
#' Beverton, R.J.H. and S.J. Holt, 1956. A review of methods for estimating mortality
#' rates in exploited fish populations, with special reference to sources of bias in
#' catch sampling. \emph{Rapports et Proces verbaux des Reunions}, Conseil Table3
#'
#' Chapman, D., and D.S Robson, 1960. The analysis of a catch curve.
#' \emph{Biometrics}, 354-368
#'
#' Edser, T., 1908. Note on the number of plaice at each length, in certain samples
#' from the southern part of the North Sea, 1906. \emph{Journal of the Royal
#' Statistical Society}, 686-690
#'
#' Heincke, F., 1913. Investigations on the plaice. General report. 1. The plaice fishery
#' and protective regulations. Part I. \emph{Rapp.P.-v.Reun.CIEM}, 17A:1-153 + Annexes
#'
#' ICES, 1981. Report of the \emph{Ad hoc} working group on the use of effort data in
#' assessment, Copenhagen, 2-6 March 1981. \emph{ICES C.M.} 1981/G:5 (mimeo)
#'
#' Jones, R., and N.P. Van Zalinge, 1981. Estimates of mortality rate and population size
#' for shrimp in Kuwait waters. \emph{Kuwait Bull. Mar. Sci}, 2, 273-288
#'
#' Pauly, D., 1983. Length-converted catch curves: a powerful tool for fisheries research
#' in the tropics (part I). \emph{ICLARM Fishbyte}, 1(2), 9-13
#'
#' Pauly, D., 1984. Length-converted catch curves: a powerful tool for fisheries
#' research in the tropics (part II). \emph{ICLARM Fishbyte}, 2(1), 17-19
#'
#' Pauly, D., 1984. Length-converted catch curves: a powerful tool for fisheries
#' research in the tropics (III: Conclusion). \emph{ICLARM Fishbyte}, 2(3), 9-10
#'
#' Ricker, W.E., 1987. Computation and interpretation of biological statistics of fish
#' populations. \emph{Bull.Fish.Res.Board Can.}, (191):382 p.
#'
#' Robson, D.S., and D.G. Chapman, 1961. Catch curves and mortality rates.
#' \emph{Trans.Am.Fish.Soc.}, 90(2):181-189
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' Van Sickle, J. 1977. Mortality rates from size distributions: the application of a
#' conservation law. \emph{Oecologia, Berl.}, 27(4):311-318
#'
#' @export

catchCurve <- function(param, catch_columns = NA, cumulative = FALSE,
                       calc_ogive = FALSE, reg_int = NULL){

  res <- param

  if("midLengths" %in% names(res)) classes <- as.character(res$midLengths)
  if("age" %in% names(res)) classes <- as.character(res$age)
  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  constant_dt <- FALSE
  if(is.na(catch_columns[1]) & (class(res$catch) == 'matrix' |
     class(res$catch) == 'data.frame')){
    writeLines("Please be aware that you provided the catch as a matrix without specifiying any columns for \n the analysis. In this case the methods applies by default the catch curve with constant \n parameter system (refer to the help file for more information).")
    flush.console()
    constant_dt <- TRUE
    catch <- res$catch
  }

  # non cumulative catch curve
  if(cumulative == FALSE){
    if(is.na(catch_columns[1]) & constant_dt == FALSE) catch <- res$catch
    if(!is.na(catch_columns[1])){
      catchmat <- res$catch[,(catch_columns)]
      if(length(catch_columns) > 1){
        catch <- rowSums(catchmat, na.rm = TRUE)
      }else catch <- catchmat
    }
  }
  #  cumulative catch curve
  if(cumulative){
    if(is.na(catch_columns[1])) catch <- rev(cumsum(rev(res$catch)))
    if(!is.na(catch_columns[1])){
      catchmat <- res$catch[,(catch_columns)]
      if(length(catch_columns) > 1){
        catchpre <- rowSums(catchmat, na.rm = TRUE)
      }else catchpre <- catchmat
      catch <- rev(cumsum(rev(catchpre)))
      #catch <- rev(cumsum(rev(res$catch[,(catch_columns)])))
    }
  }

  # Error message if catch and age do not have same length
  #   Linearised catch curve with constant time intervals
  if(constant_dt){
    if("midLengths" %in% names(res) == TRUE) stop(noquote(
      "The catch curve with constant time interval is not applicable to length-frequency data. Please provide a catch vector."))

    #if(length(classes) != length(catch[,1])) stop(noquote(
    #  "Age/length classes and catch matrix do not have the same length!"))

    if(length(classes) != length(diag(as.matrix(catch)))) writeLines("Age/length classes and the real cohort in the catch matrix \ndo not have the same length. The missing age/length \nclasses will be omitted.")

    # Aged based Catch curve
    if("age" %in% names(res) == TRUE){
      #find cohort to analyse
      real.cohort <- diag(as.matrix(catch))
      catch <- c(real.cohort, rep(NA,length(classes.num) - length(real.cohort)))
    }

  }else if(class(catch) == 'numeric'){
    if(length(classes) != length(catch)) stop(noquote(
      "Age/length classes and catch vector do not have the same length!"))
  }

  # Length converted catch curve
  if("midLengths" %in% names(res) == TRUE){
    Linf <- res$Linf
    K <- res$K
    t0 <- ifelse("t0" %in% names(res), res$t0, 0)

    if((is.null(Linf) | is.null(K))) stop(noquote(
      "You need to assign values to Linf and K for the catch curve based on length-frequency data!"))

    #calculate size class interval
    midLengths <- classes.num
    interval <- midLengths[2] - midLengths[1]

    # L and t of lower length classes
    lowerLengths <- midLengths - (interval / 2)
    t_L1 <- (t0 - (1/K)) * log(1 - (lowerLengths / Linf))

    # delta t
    dt <- rep(NA,length(midLengths))
    for(x1 in 1:(length(dt)-1)){
      dt[x1] <- t_L1[x1+1] - t_L1[x1]
    }

    # x varaible
    #ln (Linf - L)
    ln_Linf_L <- log(Linf - lowerLengths)
    # t of midlengths
    t_midL <- (t0 - (1/K)) * log(1 - (midLengths / Linf))

    # y variable
    #ln C(L1,Linf)
    lnC <- log(catch)
    # ln( Catch / delta t)
    lnC_dt <- log(catch / dt)
    lnC_dt[which(lnC_dt == -Inf)] <- NA   ### OR zero???

    if(cumulative == FALSE){
      xvar = t_midL
      yvar = lnC_dt
      xname = "t_midL"
      yname = "lnC_dt"
      xlabel = "Relative age [yrs]"
      ylabel = "ln(C/dt)"
    }

    if(cumulative){
      xvar = ln_Linf_L
      yvar = lnC
      xname = "ln_Linf_L"
      yname = "lnC"
      xlabel = "ln (Linf - L)"
      ylabel = "ln C(L, Linf)"
    }
  }

  # Aged based Catch curve
  if("age" %in% names(res) == TRUE){
    # delta t
    if(constant_dt == FALSE){
      dt <- rep(NA,length(classes.num))
      for(x1 in 1:(length(dt)-1)){
        dt[x1] <- classes.num[x1+1] - classes.num[x1]
      }
    }
    if(constant_dt) dt <- rep(classes.num[2] - classes.num[1], length(classes.num))

    # x variable
    # (t + dt) / 2   ==   x
    if(cumulative == FALSE) tplusdt_2 <- classes.num + (dt / 2)
    if(cumulative) tplusdt_2 <- classes.num

    # y variable
    #ln C(L1,Linf)
    lnC <- log(catch)
    # ln( Catch / delta t)     ==    y
    lnC_dt <- log(catch / dt)


    if(constant_dt){
      xvar = classes.num
      xname = "classes.num"
      xlabel = "Age [yrs]"
      yvar = lnC
      yname = "lnC"
      ylabel = "ln C(t, inf)"
    }

    if(cumulative == FALSE & constant_dt == FALSE){
      xvar = tplusdt_2
      xname = "tplusdt_2"
      xlabel = "Age [yrs]"
      yvar = lnC_dt
      yname = "lnC_dt"
      ylabel = "ln(C/dt)"
    }

    if(cumulative){
      xvar = tplusdt_2
      xname = "tplusdt_2"
      xlabel = "Age [yrs]"
      yvar = lnC
      yname = "lnC"
      ylabel = "ln C(t, inf)"
    }
  }

  #for plot
  #minY <- ifelse(min(yvar,na.rm=TRUE) < 0, min(yvar,na.rm=TRUE),0)
  minY <- min(yvar,na.rm=TRUE)
  maxY <- max(yvar,na.rm=TRUE) + 1
  xlims <- c(0, max(xvar,na.rm=TRUE))

  #identify plot
  if(is.null(reg_int)){
    writeLines("Please choose the minimum and maximum point in the graph \nto include for the regression line!")
    flush.console()
    dev.new()#noRStudioGD = TRUE)
    op <- par(mfrow = c(1,1),
              c(5, 4, 4, 2) + 0.1,
              oma = c(2, 1, 0, 1) + 0.1)
    plot(x = xvar,y = yvar, ylim = c(minY,maxY), xlim = xlims,
         xlab = xlabel, ylab = ylabel, type = "n")
    mtext(side = 3, "Click on two numbers. Escape to Quit.",
          xpd = NA, cex = 1.25)
    text(xvar, yvar, labels=as.character(order(xvar)), cex= 0.7)
    cutter <- identify(x = xvar, y = yvar,
                       labels = order(xvar), n=2)
    par(op)

    if(is.na(cutter[1]) | is.nan(cutter[1]) |
       is.na(cutter[2]) | is.nan(cutter[2]) ) stop(noquote("You did not choose any points in the graph. Please re-run the function and choose points in the graph!"))

    dev.off()

  }
  if(!is.null(reg_int)){
    cutter <- reg_int
  }
  if(length(cutter) != 2) stop("You have to provide 2 numbers in reg_int.")

  #calculations + model
  df.CC <- as.data.frame(cbind(xvar,yvar))
  df.CC.cut <- df.CC[cutter[1]:cutter[2],]
  lm1 <- lm(yvar ~ xvar, data = df.CC.cut)
  sum_lm1 <- summary(lm1)
  r_lm1 <- sum_lm1$r.squared
  intercept_lm1 <- sum_lm1$coefficients[1]
  slope_lm1 <- sum_lm1$coefficients[2]
  se_slope_lm1 <- sum_lm1$coefficients[4]


  #fit of regression line
  lm1.fit <- sum_lm1$r.squared
  Z_lm1 <- abs(slope_lm1)
  SE_Z_lm1 <- abs(se_slope_lm1)
  confi <-  abs(se_slope_lm1) * qt(0.975,sum_lm1$df[2])
  conf_Z_lm1 <- Z_lm1 + c(-confi,confi)


  # special case when cumulative and length-frequency data
  if(cumulative & "midLengths" %in% names(res) == TRUE){
    Z_lm1 <- Z_lm1 * K
    SE_Z_lm1 <- SE_Z_lm1 * K
  }

  #save all in list
  ret <- c(res,list(
    xvar = xvar,
    yvar = yvar,
    reg_int = cutter,
    linear_mod = lm1,
    Z =  Z_lm1,
    se = SE_Z_lm1,
    confidenceInt = conf_Z_lm1
  ))
  names(ret)[names(ret) == "xvar"] <- xname
  names(ret)[names(ret) == "yvar"] <- yname
  class(ret) <- "catchCurve"

  # Calculate selection ogive from catch curve and add to ret
  if(calc_ogive & cumulative) stop(noquote("It is not possible to estimate the selection ogive for the cumulative catch curve."))
  if(calc_ogive){

    # only use part of catch and t which is not fully exploited by the gear
    t_ogive <- xvar[1:(cutter[1]-1)]
    dt_ogive <- dt[1:(cutter[1]-1)]
    if("age" %in% names(res) == TRUE &
       class(catch) == 'matrix' | class(catch) == 'data.frame'){
      catch_ogive <- catch[1:(cutter[1]-1)] ## catch.cohort
    }else catch_ogive <- catch[1:(cutter[1]-1)]

    # calculate observed selection ogive
    Sobs <- catch_ogive/(dt_ogive * exp(intercept_lm1 - Z_lm1*t_ogive))

    # dependent vairable in following regression analysis
    ln_1_S_1 <- log((1/Sobs) - 1)

    # get rid of Inf
    ln_1_S_1[which(ln_1_S_1 == Inf)] <- NA
    t_ogive[which(t_ogive == Inf)] <- NA

    #regression analysis to caluclate T1 and T2
    mod_ogive <- lm(ln_1_S_1 ~ t_ogive, na.action = na.omit)
    sum_lm_ogive <- summary(mod_ogive)
    T1 <- sum_lm_ogive$coefficients[1]
    T2 <- abs(sum_lm_ogive$coefficients[2])

    # calculate estimated selection ogive
    Sest <- 1/(1+exp(T1 - T2*xvar))

    # selection parameters
    t50 <- T1/T2
    t75 <- (T1 + log(3))/T2
    if(!is.null(res$Linf) & !is.null(res$K)){
      if(is.null(res$t0)) t0 = 0
      L50 <- Linf*(1-exp(-K*(t50-t0)))
      L75 <- Linf*(1-exp(-K*(t75-t0)))
    }

    ret2 <- c(ret,list(
      intercept = intercept_lm1,
      linear_mod_sel = mod_ogive,
      Sobs = Sobs,
      ln_1_S_1 = ln_1_S_1,
      Sest = Sest,
      t50 = t50,
      t75 = t75))
    if(exists("L50")) ret2$L50 = L50
    if(exists("L75")) ret2$L75 = L75
    if(exists("L50")) names(ret2)[which(ret2 %in% L50)] <- "L50"
    if(exists("L75")) names(ret2)[which(ret2 %in% L75)] <- "L75"

    class(ret2) <- "catchCurve"
    plot(ret2, plot_selec=TRUE)
    return(ret2)
  }else {plot(ret)
    return(ret)}
}
