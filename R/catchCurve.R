#' @title Catch curve
#'
#' @description  The funciton estimates the total mortality (Z) from age composition
#'    or length-frequency data with the linearized catch curve or the linearized
#'    length converted catch curve, respectively. It allows for applying the
#'    cumulative catch curve and to estimate gear selectivity.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \code{midLengths} or \code{age}: midpoints of the length class as vector (length-frequency
#'   data) or ages as vector (age composition data),
#'   \item \code{Linf}: infinite length for investigated species in cm [cm],
#'   \item \code{K}: growth coefficent for investigated species per year [1/year],
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \code{catch}: catch as vector, or a matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#' }
#' @param catch_column numerical; indicating the column of the catch matrix which should be
#'   used for the analysis.
#' @param cumulative logical; if \code{TRUE} instead of normal catch curve the cumulative
#'   catch curve is applied (Jones and van Zalinge method)
#' @param calc_ogive logical; if \code{TRUE} the selection ogive is additionally
#'   calculated from the catch curve
#'
#' @keywords function mortality Z catchCurve
#'
#' @examples
#' \donttest{
#' #_______________________________________________
#' # Variable paramter system (with catch vector)
#' # based on length frequency data
#' # load data
#' data(goatfish)
#'
#' # run model
#' catchCurve(goatfish)
#'
#' # based on age composition data
#' # load data
#' data(whiting)
#'
#' # run model
#' catchCurve(param=whiting, catch_column = 1)
#'
#' #_______________________________________________
#' # Constant parameter system based on age composition data (with catch matrix)
#' catchCurve(param = whiting)
#'
#'
#' #_______________________________________________
#' # Cumulative Catch Curve
#' # based on length frequency data
#' # load data
#' data(goatfish)
#'
#' # run model
#' catchCurve(param = goatfish, cumulative = TRUE)
#'
#'
#' # based on age composition data
#' data(synCAA2)
#'
#' # run model
#' catchCurve(synCAA2, cumulative = TRUE)
#'
#'
#' #_______________________________________________
#' # Catch Curve with estimation of selection ogive
#' # load data
#' data(synLFQ3)
#'
#' # run model
#' catchCurve(synLFQ3, calc_ogive = TRUE)
#'
#'  }
#'
#' @details This function includes a so called 'locator' function, which asks you to
#'   choose two points from a graph manually. The two points which you choose by clicking
#'   on the plot in the graphical device represent the start and end of the data points,
#'   which should be used for the analysis. Based on these points the regression line
#'   is calculated.
#'   When the selection ogive
#'   is calculated by means of the catch curve the assumption is made, that Z is constant
#'   for all year classes or length groups respectively. Accoring to Sparre and Venema
#'   (1998) this assumption might be true, because F is smaller for young fish
#'   (Selectivity) while M is higher for young fish (high natural mortality). The selectivity
#'   for not fully exploited old fish (e.g. due to gillnet fishery) can not be calculated yet
#'   by use of the catch curve.
#'   Based on the format of the list argument "catch" the function automatically
#'   distinguishes between the catch curve with variable parameter system (if catch is a
#'   vector) and the one with constant parameter system (if catch is a matrix or a
#'   data.frame). In the case of the variable parameter system the catches of one year are
#'   assumed to represent the catches during the entire life span of a so called
#'   pseudo-cohort.
#'   The cumulative catch curve does not allow for the estimation of the selectivity
#'   ogive.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{tplusdt_2} or \strong{t_midL}: relative ages depending on input
#'      vector (\code{age} or \code{midLengths}),
#'   \item \strong{lnC_dt}: rearranged catches,
#'   \item \strong{reg_int}: the interval used for the regression analysis,
#'   \item \strong{Z}: the total mortality,
#'   \item \strong{se}: the standard error of the total mortality;}
#' in case of calc_ogive, additionally:
#' \itemize{
#'   \item \strong{intercept}: intercept of regression analysis,
#'   \item \strong{Sobs}: observed selection curve,
#'   \item \strong{ln_1_S_1}: dependent variable of regression analysis for
#'   selectivity parameters,
#'   \item \strong{Sest}: estimated selection curve,
#'   \item \strong{t50}: age at first capture (age at which fish have a 50%
#'   probability to be caught),
#'   \item \strong{t75}: age at which fish have a 75% probability to be caught,
#'   \item \strong{L50}: length at first capture (length at which fish have a 50%
#'   probability to be caught),
#'   \item \strong{L75}: length at which fish have a 75% probability to be caught;
#' }
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

catchCurve <- function(param, catch_column = NA, cumulative = FALSE,
                       calc_ogive = FALSE){

  res <- param

  #IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#
  #      NON CUMULATIVE CATCH CURVE   #
  #IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#
  if(cumulative == FALSE){
    if(is.na(catch_column)) catch <- res$catch
    if(!is.na(catch_column)) catch <- res$catch[,catch_column]
    if("midLengths" %in% names(res)) classes <- as.character(res$midLengths)
    if("age" %in% names(res)) classes <- as.character(res$age)

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
    }

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])


    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #   Linearised catch curve with constant time intervals    #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){

      if("midLengths" %in% names(res) == TRUE) stop("The catch curve with constant time interval is not applicable to length frequency data.")

      #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
      #   Aged based Catch curve   #
      #HHHHHHHHHHHHHHHHHHHHHHHHHHHH#
      if("age" %in% names(res) == TRUE){

        #find cohort to analyse
        real.cohort <- diag(as.matrix(catch))
        catch.cohort <- c(real.cohort,
                          rep(NA,length(classes.num) - length(real.cohort)))

        # ln( Catch )     ==    y
        lnC <- log(catch.cohort)

        #for plot
        minlnC <- ifelse(min(lnC,na.rm=TRUE) < 0, min(lnC,na.rm=TRUE),0)
        maxlnC <- max(lnC,na.rm=TRUE) + 1

        #identify plot
        op <- par(mfrow=c(1,1))
        plot(x = classes.num,y = lnC, ylim = c(minlnC,maxlnC),
             xlab = "Age [yrs]", ylab = "ln(C)")
        print("Please choose the minimum and maximum point in the graph to include for the regression line!")
        cutter <- identify(x = classes.num, y = lnC,
                           labels = order(classes.num), n=2)
        par(op)

        #calculations + model
        df.CC <- as.data.frame(cbind(classes.num,lnC))
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


        #save all in list
        ret <- c(res,list(
          classes.num = classes.num,
          lnC = lnC,
          reg_int = cutter,
          Z =  Z_lm1,
          se = SE_Z_lm1
        ))
        class(ret) <- "catchCurve"

        # Calculate selection ogive from catch curve and add to ret
        if(calc_ogive){

          dt <- rep(classes.num[2]-classes.num[1],length(classes.num))

          # only use part of catch and t which is not fully exploited by the gear
          t_ogive <- classes.num[1:(cutter[1]-1)]   # replace t_midL by universal object
          dt_ogive <- dt[1:(cutter[1]-1)]
          catch_ogive <- catch.cohort[1:(cutter[1]-1)]

          # calculate observed selection ogive
          Sobs <- catch_ogive/(dt_ogive * exp(intercept_lm1 - Z_lm1*t_ogive))

          # dependent vairable in following regression analysis
          ln_1_S_1 <- log((1/Sobs) - 1)

          #regression analysis to caluclate T1 and T2
          sum_lm_ogive <- summary(lm(ln_1_S_1 ~ t_ogive, na.action = na.omit))
          T1 <- sum_lm_ogive$coefficients[1]
          T2 <- abs(sum_lm_ogive$coefficients[2])

          # calculate estimated selection ogive
          Sest <- 1/(1+exp(T1 - T2*classes.num))

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
          #plot(ret2,plot.selec=TRUE)
          return(ret2)
        }else plot(ret) ; return(ret)
      }
    }

    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    #    Linearised catch curve with variable time intervals   #
    #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if(class(catch) == 'numeric'){

      #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
      #      Length converted catch curve     #
      #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
      if("midLengths" %in% names(res) == TRUE){

        Linf <- res$Linf
        K <- res$K
        t0 <- ifelse("t0" %in% names(res), res$t0, 0)

        if((is.null(Linf) | is.null(K))) stop("You need to assign values to Linf and K for the Catch curve based on length frequency data!")

        #calculate size class interval
        midLengths <- classes.num
        interval <- midLengths[2] - midLengths[1]

        # t of lower length classes
        lowerLengths <- midLengths - (interval / 2)
        t_L1 <- (t0 - (1/K)) * log(1 - (lowerLengths / Linf))

        # delta t
        dt <- rep(NA,length(midLengths))
        for(x1 in 1:(length(dt)-1)){
          dt[x1] <- t_L1[x1+1] - t_L1[x1]
        }

        # t of midlengths
        t_midL <- (t0 - (1/K)) * log(1 - (midLengths / Linf))
        # ln( Catch / delta t)
        lnC_dt <- log(catch / dt)

        lnC_dt[which(lnC_dt == -Inf)] <- NA   ### OR zero???

        #for plot
        minlnC_dt <- ifelse(min(lnC_dt,na.rm=TRUE) < 0, min(lnC_dt,na.rm=TRUE),0)
        maxlnC_dt <- max(lnC_dt,na.rm=TRUE) + 1

        #identify plot
        op <- par(mfrow=c(1,1))
        plot(x = t_midL,y = lnC_dt, ylim = c(minlnC_dt,maxlnC_dt),
             xlab = "Relative age [yrs]", ylab = "ln(C/dt)")
        print("Please choose the minimum and maximum point in the graph to include for the regression line!")
        cutter <- identify(x = t_midL, y = lnC_dt,
                           labels = order(t_midL), n=2)
        par(op)

        #calculations + model
        df.CC <- as.data.frame(cbind(t_midL,lnC_dt))
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

        #save all in list
        ret <- c(res,list(
          t_midL = t_midL,
          lnC_dt = lnC_dt,
          reg_int = cutter,
          Z =  Z_lm1,
          se = SE_Z_lm1
        ))
        class(ret) <- "catchCurve"

        # Calculate selection ogive from catch curve and add to ret
        if(calc_ogive){

          # only use part of catch and t which is not fully exploited by the gear
          t_ogive <- t_midL[1:(cutter[1]-1)]   # replace t_midL by universal object
          dt_ogive <- dt[1:(cutter[1]-1)]
          catch_ogive <- catch[1:(cutter[1]-1)]

          # calculate observed selection ogive
          Sobs <- catch_ogive/(dt_ogive * exp(intercept_lm1 - Z_lm1*t_ogive))

          # dependent vairable in following regression analysis
          ln_1_S_1 <- log((1/Sobs) - 1)

          #regression analysis to caluclate T1 and T2
          sum_lm_ogive <- summary(lm(ln_1_S_1 ~ t_ogive))
          T1 <- sum_lm_ogive$coefficients[1]
          T2 <- abs(sum_lm_ogive$coefficients[2])

          # calculate estimated selection ogive
          Sest <- 1/(1+exp(T1 - T2*t_midL))

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
            Sobs = Sobs,
            ln_1_S_1 = ln_1_S_1,
            Sest = Sest,
            t50 = t50,
            t75 = t75))
          if(exists("L50")){
            ret2$L50 = L50
            names(ret2)[which(ret2 %in% L50)] <- "L50"
          }
          if(exists("L75")){
            ret2$L75 = L75
            names(ret2)[which(ret2 %in% L75)] <- "L75"
          }

          class(ret2) <- "catchCurve"
          plot(ret2, plot.selec=TRUE)
          return(ret2)
        }else plot(ret) ; return(ret)
      }

      #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
      #    Aged based Catch curve    #
      #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
      if("age" %in% names(res) == TRUE){

        # delta t
        dt <- rep(NA,length(classes.num))
        for(x1 in 1:(length(dt)-1)){
          dt[x1] <- classes.num[x1+1] - classes.num[x1]
        }

        # (t + dt) / 2   ==   x
        tplusdt_2 <- classes.num + (dt / 2)

        # ln( Catch / delta t)     ==    y
        lnC_dt <- log(catch / dt)


        #for plot
        minlnC_dt <- ifelse(min(lnC_dt,na.rm=TRUE) < 0, min(lnC_dt,na.rm=TRUE),0)
        maxlnC_dt <- max(lnC_dt,na.rm=TRUE) + 1

        #identify plot
        op <- par(mfrow=c(1,1))
        plot(x = tplusdt_2,y = lnC_dt, ylim = c(minlnC_dt,maxlnC_dt),
             xlab = "Age [yrs]", ylab = "ln(C/dt)")
        print("Please choose the minimum and maximum point in the graph to include for the regression line!")
        cutter <- identify(x = tplusdt_2, y = lnC_dt,
                           labels = order(tplusdt_2), n=2)
        par(op)

        #calculations + model
        df.CC <- as.data.frame(cbind(tplusdt_2,lnC_dt))
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

        #save all in list
        ret <- c(res,list(
          tplusdt_2 = tplusdt_2,
          lnC_dt = lnC_dt,
          reg_int = cutter,
          Z =  Z_lm1,
          se = SE_Z_lm1
        ))
        class(ret) <- "catchCurve"

        # Calculate selection ogive from catch curve and add to ret
        if(calc_ogive){

          # only use part of catch and t which is not fully exploited by the gear
          t_ogive <- t_midL[1:(cutter[1]-1)]   # replace t_midL by universal object
          dt_ogive <- dt[1:(cutter[1]-1)]
          catch_ogive <- catch[1:(cutter[1]-1)]

          # calculate observed selection ogive
          Sobs <- catch_ogive/(dt_ogive * exp(intercept_lm1 - Z_lm1*t_ogive))

          # dependent vairable in following regression analysis
          ln_1_S_1 <- log((1/Sobs) - 1)

          #regression analysis to caluclate T1 and T2
          sum_lm_ogive <- summary(lm(ln_1_S_1 ~ t_ogive))
          T1 <- sum_lm_ogive$coefficients[1]
          T2 <- abs(sum_lm_ogive$coefficients[2])

          # calculate estimated selection ogive
          Sest <- 1/(1+exp(T1 - T2*tplusdt_2))

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
            Sobs = Sobs,
            ln_1_S_1 = ln_1_S_1,
            Sest = Sest,
            t50 = t50,
            t75 = t75))
          if(exists("L50")){
            ret2$L50 = L50
            names(ret2)[which(ret2 %in% L50)] <- "L50"
          }
          if(exists("L75")){
            ret2$L75 = L75
            names(ret2)[which(ret2 %in% L75)] <- "L75"
          }

          class(ret2) <- "catchCurve"
          plot(ret2,plot.selec=TRUE)
          return(ret2)
        }else plot(ret) ; return(ret)

        # plot results
        plot(ret)

        return(ret)
      }
    }
  }


  #IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#
  #       CUMULATIVE CATCH CURVE      #
  #IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII#
  if(cumulative){
    if(is.na(catch_column)) cumCatch <- rev(cumsum(rev(res$catch)))
    if(!is.na(catch_column)) cumCatch <- rev(cumsum(rev(res$catch[,catch_column])))

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

        if(is.null(Linf) | is.null(K)) stop("You need to assign values to Linf and K for the cumulated catch curve based on length frequency data!")

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
        minlnC <- ifelse(min(ln_C,na.rm=TRUE) < 0, min(ln_C,na.rm=TRUE),0)
        maxlnC <- max(ln_C,na.rm=TRUE) + 1

        #identify plot
        writeLines("Please choose the minimum and maximum point in the graph \nto include for the regression line! Then press 'Finish'!")
        if(.Platform$OS.type == "unix") quartz()
        if(.Platform$OS.type == "windows") windows()
        op <- par(mfrow = c(1,1),
                  c(5, 4, 4, 2) + 0.1)
        plot(x = ln_Linf_L,y = ln_C, ylim = c(minlnC,maxlnC),
             xlab = "ln (Linf - L)", ylab = "ln C(L, Linf)")
        cutter <- identify(x = ln_Linf_L, y = ln_C,
                           labels = order(ln_Linf_L), n=2)
        par(op)

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

        #save all in list
        ret <- c(res,list(
          classes.num = classes.num,
          ln_Linf_L = ln_Linf_L,
          ln_C = ln_C,
          reg_int = cutter,
          Z =  Z_lm1,
          se = SE_Z_lm1
        ))

        class(ret) <- "catchCurve"

        # plot results
        plot(ret)

        return(ret)
      }

      #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
      #    Aged based Catch curve    #
      #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
      if("age" %in% names(res) == TRUE){

        tplusdt_2 <- classes.num

        #ln C(L1,Linf)
        ln_C <- log(cumCatch)

        #for plot
        minlnC <- ifelse(min(ln_C,na.rm=TRUE) < 0, min(ln_C,na.rm=TRUE),0)
        maxlnC <- max(ln_C,na.rm=TRUE) + 1

        #identify plot
        writeLines("Please choose the minimum and maximum point in the graph \nto include for the regression line! Then press 'Finish'!")
        if(.Platform$OS.type == "unix") quartz()
        if(.Platform$OS.type == "windows") windows()
        op <- par(mfrow = c(1,1),
                  c(5, 4, 4, 2) + 0.1)
        plot(x = tplusdt_2,y = ln_C, ylim = c(minlnC,maxlnC),
             xlab = "Age [yrs]", ylab = "ln C(t, inf)")
        cutter <- identify(x = tplusdt_2, y = ln_C,
                           labels = order(tplusdt_2),n=2)
        par(op)

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

        #save all in list
        ret <- c(res,list(
          tplusdt_2 = tplusdt_2,
          ln_C = ln_C,
          reg_int = cutter,
          Z =  Z_lm1,
          se = SE_Z_lm1
        ))
        class(ret) <- "catchCurve"

        # plot results
        plot(ret)

        return(ret)
      }
    }
  }
}
