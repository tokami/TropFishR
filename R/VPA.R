#' @title Virtual Population Analysis (VPA) model
#'
#' @description This function applies the Virtual Population Analysis or Cohort analysis,
#'    respectively. A method used to estimate fishing mortality per age/length group and
#'    to get a first estimate of stock biomass.
#'
#' @param param A list consisting of following parameters:
#'   \code{$age} or \code{$midLengths} midpoints of the length class as vector (length frequency
#'   data) or ages as vector (age composition data),
#'   \code{$Linf} Infinite length for investigated species in cm [cm],
#'   \code{$K} Growth coefficent for investigated species per year [1/year],
#'   \code{t0} Theoretical time zero, at which individuals of this species hatch,
#'   \code{M} Natural mortality [1/year],
#'   \code{a} length-weight relationship coefficent (W = a * L^b),
#'   \code{b} length-weight relationship coefficent (W = a * L^b),
#'   \code{catch} Catch as vector, or a matrix with catches of subsequent years if
#'   the catch curve with constat time intervals should be applied;
#' @param terminalF A fishing mortality value which is used as the terminal FM for the last age/length group.
#' @param analysis.type Determines which type of assessment should be done,
#'   options: "VPA" for classical age-based VPA, "CA" for age- or length-based
#'   Cohort analysis
#' @param catch_corFac optional: Correction factor for catch, in case provided
#'   catch does spatially or temporarily not reflect catch for fishing ground of
#'   a whole year.
#' @param algorithm An Algorithm to use to solve for fishing mortality. The default
#'   setting \code{algorithm="new"} uses \code{\link[stats]{optimize}},
#'   while \code{algorithm="old"} uses the algorithm described by
#'   Sparre and Venema (1998).
#'
#' @details VPA and Cohort analyses
#'
#' @examples
#' # Virtual Popuation Analysis with age-composition data
#' data(whiting)
#' output <- VPA(whiting, terminalF = 0.5, analysis.type = "VPA")
#'
#' # Pope's Cohort Analysis with age-composition data
#' data(whiting)
#' output <- VPA(whiting, terminalF = 0.5, analysis.type = "CA")
#'
#' # Jones's Cohort Analysis with length-composition data
#' data(hake)
#' output <- VPA(hake, terminalF = 0.5, analysis.type = "CA")
#'
#' @references
#' Jones?
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' References for weight-length relationship parameters (a & b):
#' Dorel, D., 1986. Poissons del'Atlantique nord-est relations taille-poids. Institut Francais de Recherche
#' pour l'Exploitation de la Mer. Nantes, France. 165 p.
#'
#' @export

VPA <- function(param, terminalF, analysis.type, catch_corFac = NA, algorithm="new"){

  res <- param
  catch <- res$catch

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #           AGE BASED VPA AND COHORT ANALYSIS              #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#

  if("age" %in% names(res) == TRUE){

    classes <- as.character(res$age)

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Ages and catch do not have the same length!")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Ages and catch do not have the same length!")
    }

    a <- res$a
    b <- res$b
    M <- res$M

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      #find cohort to analyse
      real.cohort <- diag(as.matrix(catch))      ##CHECK!!! TAKES ALWAYS THE FIRST OBSERVATION IN FIRST COLUMN TO START FINDING REAL COHORT
      catch.cohort <- c(real.cohort,
                        rep(NA,length(classes.num) - length(real.cohort)))
    }
    if(class(catch) == 'numeric'){
      catch.cohort <- catch
    }

    #Correct catch  if not representative for one year
    if(!is.na(catch_corFac)) catch_cor <- catch.cohort * catch_corFac
    if(is.na(catch_corFac)) catch_cor <- catch.cohort

    #Survivors    #N(L1)=(N(L2)*H(L1,L2)+C(L1,L2)) *H(L1,L2)
    survivors <- rep(NA,length(classes.num))

    # survivors last size class
    lastLengthClass <- max(which(!is.na(catch_cor)),na.rm=TRUE)  ### CHECK!!!!: EXCLUDES THE LAST CLASSES WHICH ARE NA IN CATCH DUE TO TOO LESS YEARS SAMPLED!
    survivors[lastLengthClass] <-
      catch_cor[lastLengthClass] / ((terminalF/(terminalF + M)) * (1 - exp(-(terminalF + M))))


    #   Age-based Cohort Analysis (Pope's cohort analysis)
    if(analysis.type == "CA"){
      # other survivors
      for(x3 in (lastLengthClass-1):1){
        survivors[x3] <- (survivors[x3+1] * exp((M/2)) +
                            catch_cor[x3] ) * exp((M/2))
      }

      #F
      FM <- rep(NA,length(classes.num))
      FM[lastLengthClass] <- terminalF
      for(x5 in 1:(lastLengthClass-1)){
        FM[x5] <- log(survivors[x5]/survivors[x5+1]) - M
      }
    }

    # Traditional VPA
    if(analysis.type == "VPA"){
      #other survivors and fishing mortality
      ###IMPROVABLE BY MAKING THE STEP CHOOSABLE, MEANING THE USER CAN CHOOSE THE RESOLUTION
      FM <- rep(NA,length(classes.num))
      FM[lastLengthClass] <- terminalF

      for(num_class in (lastLengthClass-1):1){

        sur.C <- catch_cor[num_class]
        sur.Ntplus1 <- survivors[(num_class+1)]
        sur.M <- M
        LHS <-  sur.C / sur.Ntplus1
        sur.F <- 0
        seqi <- c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7)

        if(algorithm == "old"){
          LHS <-  sur.C / sur.Ntplus1
          sur.F <- 0
          seqi <- c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7)
          #trail and error
          for(y in seqi){
            stepi <- y
            for(x in seq(sur.F,10,stepi)){
              sur.F <- x
              RHS <- (sur.F/(sur.F + sur.M)) * (exp(sur.F+sur.M) - 1)
              if(LHS-RHS < 0) break
            }
            sur.F = x-stepi
          }
        }

        if(algorithm == "new"){
          Fcalc <- function(sur.F=sur.M){
            ((sur.F/(sur.F+sur.M)) * (exp(sur.F+sur.M) - 1) - (sur.C / sur.Ntplus1))^2
          }
          tmp <- optimize(Fcalc, interval=c(0,100))
          sur.F <- tmp$min
        }

        #fill F
        FM[num_class] <- sur.F

        #fill survivors
        survivors[num_class] <- survivors[(num_class+1)] *
          exp(sur.F + sur.M)
      }
    }

    # Z
    Z <- rep(NA,length(classes.num))
    for(x6 in 1:(length(Z))){
      Z[x6] <- M  +  FM[x6]
    }

    #Annual mean Nr
    annualMeanNr <- rep(NA,length(classes.num))
    for(x7 in 1:(length(annualMeanNr-1))){
      annualMeanNr[x7] <- (survivors[x7] -
                             survivors[x7+1]) / Z[x7]
    }

    #Mean body weight
    meanBodyWeight <- a * classes.num ^ b

    #Mean biomass
    meanBiomass <- annualMeanNr * meanBodyWeight
    meanBiomassTon <- meanBiomass/1000

    #Yield
    yield <- catch_cor * meanBodyWeight
    yieldTon <- yield/1000

    #FOR PLOT
    #Survivors rearranged
    survivors_rea <- rep(NA,length(classes.num))
    for(x8 in 1:(length(survivors_rea)-1)){
      survivors_rea[x8] <- survivors[x8+1]
    }
    survivors_rea[length(survivors_rea)] <- 0


    #Calculate natural losses
    natLoss <- rep(NA,length(classes.num))
    for(x9 in 1:length(natLoss)){
      natLoss[x9] <- survivors[x9] - survivors_rea[x9] -
        catch_cor[x9]
    }

    #put together in dataframe
    df.VPAnew <- data.frame(survivors = survivors_rea,
                            nat.losses = natLoss,
                            catch = catch_cor)

    #transpose matrix for barplot function
    df.VPAnew <- t(as.matrix(df.VPAnew))
    colnames(df.VPAnew) <- classes.num

    #save x axis positions
    max_sur <- round(max(survivors,na.rm=TRUE),digits=0)
    dim_sur <- 10 ^ (nchar(max_sur)-1)
    max_FM <- ceiling(max(FM,na.rm=TRUE))
    max_clas <- max(classes.num,na.rm=TRUE)
    par(new = FALSE)
    mids <- barplot(df.VPAnew, xlab="", ann=TRUE,
                    ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur))

    #create VPA plot
    par(mar = c(5, 4, 4, 4) + 0.3)
    barplot(df.VPAnew,col=c('darkgreen','purple','yellow'),
            xlab = "Age", ylab = "Population",xlim=c(0,ceiling(max(mids))),
            ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur))
    legend(x=mids[(which(classes.num == max_clas)-2)],
           y = ceiling(max_sur/dim_sur)*dim_sur,
           legend = c(rownames(df.VPAnew),"Fishing mortality"),
           col = c('darkgreen','purple','yellow','red'),xpd = TRUE,
           pch=c(rep(15,3),NA), lty = c(NA,NA,NA,1), lwd=2,seg.len = 0.3,
           pt.cex = 2, x.intersp = c(0.3,0.3,0.3,0.3),merge=TRUE,
           y.intersp = 0.6, box.lty=0,cex=0.8,xjust = 0,yjust = 0.8)
    par(new = TRUE,mar = c(5, 4, 4, 4) + 0.3)
    plot(mids, FM, col='red',xlim=c(0,ceiling(max(mids))),
         type = "n",axes = FALSE, bty = "n", xlab = "", ylab = "",ann=TRUE)
    lines(x=mids,y=FM,col='red',lwd=2)
    usr <- par("usr")
    par(usr=c(usr[1:2], 0, max_FM))
    axis(4,at=pretty(c(0,max_FM)))
    mtext("fishing mortatlity", side=4, line=3)

    #save all in list
    ret <- c(res,list(
      classes.num = classes.num,
      catch.cohort = catch.cohort,
      FM = FM,
      Z = Z,
      survivors = survivors,
      annualMeanNr = annualMeanNr,
      meanBodyWeight = meanBodyWeight,
      meanBiomassTon = meanBiomassTon,
      yieldTon = yieldTon,
      natLoss = natLoss,
      plot_mat = df.VPAnew))

    class(ret) <- "VPA"
    return(ret)
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #           Jones' Length-based Cohort Analysis            #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  if("midLengths" %in% names(res) == TRUE &
     (class(catch) == 'matrix' |
      class(catch) == 'data.frame')) stop("The length-based Cohort analysis
                                          is not applicable to length frequency data.
                                          Please provide catch as vector.")

  if("midLengths" %in% names(res) == TRUE & analysis.type == "VPA") stop("Please choose
                                                         analysis.type = 'CA' for
                                                         length composition data!")

  if(class(catch) == 'numeric' & "midLengths" %in% names(res) == TRUE &
     analysis.type == "CA"){

    classes <- as.character(res$midLengths)

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Midlengths and catch do not have the same length!")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Midlengths and catch do not have the same length!")
    }

    Linf <- res$Linf
    K <- res$K
    t0 <- ifelse(is.null(res$t0),0,res$t0)
    a <- res$a
    b <- res$b
    M <- res$M

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    #calculate size class interval
    interval <- classes.num[2] - classes.num[1]

    # t of lower length classes
    lowerLength <- classes.num - (interval / 2)
    if(!is.na(catch_corFac)) catch_cor <- catch * catch_corFac
    if(is.na(catch_corFac)) catch_cor <- catch
    t_L1 <- (t0 - (1/K)) * log(1 - (lowerLength / Linf))

    # delta t
    dt <- rep(NA,length(classes.num))
    for(x1 in 1:(length(dt)-1)){
      dt[x1] <- t_L1[x1+1] - t_L1[x1]
    }

    # t of midlengths
    t_midL <- (t0 - (1/K)) * log(1 - (classes.num / Linf))

    # H (L1,L2)   #H(L1,L2)=((Linf-L1)/Linf-L2)^(M/2K)
    H <- rep(NA,length(classes.num))
    for(x2 in 1:(length(H)-1)){
      H[x2] <- ((Linf - lowerLength[x2]) /
                         (Linf - lowerLength[x2+1])) ^
        (M / (2*K))
    }

    #Survivors    #N(L1)=(N(L2)*H(L1,L2)+C(L1,L2)) *H(L1,L2)
    survivors <- rep(NA,length(classes.num))

    # survivors last size class
    survivors[length(survivors)] <-
      catch_cor[length(survivors)] / (terminalF/(terminalF + M))
    # other survivors
    for(x3 in (length(survivors)-1):1){
      survivors[x3] <- (survivors[x3+1] *
                                 H[x3] + catch_cor[x3] ) *
        H[x3]
    }

    # F/Z  #F(L1,L2)/Z(L1,L2)=C(L1,L2)/(N(L1)-N(L2))
    F_Z <- rep(NA,length(classes.num))
    for(x4 in 1:(length(F_Z)-1)){
      F_Z[x4] <- catch_cor[x4] /
        (survivors[x4] - survivors[x4+1])
    }
    F_Z[length(F_Z)] <- terminalF / (terminalF + M)

    #F  # F = M * (F_Z / 1-F_Z)
    FM <- rep(NA,length(classes.num))
    for(x5 in 1:(length(FM))){
      FM[x5] <- M  *  (F_Z[x5] / (1 - F_Z[x5]))
    }

    # Z
    Z <- rep(NA,length(classes.num))
    for(x6 in 1:(length(Z))){
      Z[x6] <- M  +  FM[x6]
    }

    #Annual mean Nr
    annualMeanNr <- rep(NA,length(classes.num))
    for(x7 in 1:(length(annualMeanNr-1))){
      annualMeanNr[x7] <- (survivors[x7] -
                                    survivors[x7+1]) / Z[x7]
    }

    #Mean body weight
    meanBodyWeight <- a * classes.num ^ b

    #Mean biomass
    meanBiomass <- annualMeanNr * meanBodyWeight
    meanBiomassTon <- meanBiomass/1000

    #Yield
    yield <- catch_cor * meanBodyWeight
    yieldTon <- yield/1000

    #FOR PLOT
    #Survivors rearranged
    survivors_rea <- rep(NA,length(classes.num))
    for(x8 in 1:(length(survivors_rea)-1)){
      survivors_rea[x8] <- survivors[x8+1]
    }
    survivors_rea[length(survivors_rea)] <- 0

    #Calculate natural losses
    natLoss <- rep(NA,length(classes.num))
    for(x9 in 1:length(natLoss)){
      natLoss[x9] <- survivors[x9] - survivors_rea[x9] -
        catch_cor[x9]
    }

    #put together in dataframe
    df.VPAnew <- data.frame(survivors = survivors_rea,
                            nat.losses = natLoss,
                            catch = catch_cor)

    #transpose matrix for barplot function
    df.VPAnew <- t(as.matrix(df.VPAnew))
    colnames(df.VPAnew) <- classes.num

    #save x axis positions
    max_sur <- round(max(survivors,na.rm=TRUE),digits=0)
    dim_sur <- 10 ^ (nchar(max_sur)-1)
    max_FM <- ceiling(max(FM,na.rm=TRUE))
    max_clas <- max(classes.num,na.rm=TRUE)
    par(new = FALSE)
    mids <- barplot(df.VPAnew, xlab="",
                    ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur))

    #create VPA plot
    par(mar = c(5, 4, 4, 4) + 0.3)
    barplot(df.VPAnew,col=c('darkgreen','purple','yellow'),
            xlab = "Midlength [cm]", ylab = "Population",xlim=c(0,ceiling(max(mids))),
            ylim = c(0,ceiling(max_sur/dim_sur)*dim_sur))
    legend(x=mids[(which(classes.num == max_clas)-3)],
           y = (ceiling(max_sur/dim_sur)*dim_sur),
           legend = c(rownames(df.VPAnew),"Fishing mortality"),
           col = c('darkgreen','purple','yellow','red'),xpd = TRUE,
           pch=c(rep(15,3),NA), lty = c(NA,NA,NA,1), lwd=2,seg.len = 0.3,
           pt.cex = 2, x.intersp = c(0.3,0.3,0.3,0.3),merge=TRUE,
           y.intersp = 0.6, box.lty=0,cex=0.8,xjust = 0,yjust = 0.8)
    par(new = TRUE,mar = c(5, 4, 4, 4) + 0.3)
    plot(mids, FM, col='red',xlim=c(0,ceiling(max(mids))),
         type = "n",axes = FALSE, bty = "n", xlab = "", ylab = "",ann=TRUE)
    lines(x=mids,y=FM,col='red',lwd=2)
    usr <- par("usr")
    par(usr=c(usr[1:2], 0, max_FM))
    axis(4,at=pretty(c(0,max_FM)))
    mtext("fishing mortatlity", side=4, line=3)
    plot1 <- recordPlot()

    #save all in list
    ret <- c(res,list(
      classes.num = classes.num,
      FM = FM,
      Z = Z,
      survivors = survivors,
      annualMeanNr = annualMeanNr,
      meanBodyWeight = meanBodyWeight,
      meanBiomassTon = meanBiomassTon,
      yieldTon = yieldTon,
      natLoss = natLoss,
      plot_mat = df.VPAnew))

    class(ret) <- "VPA"
    return(ret)
  }
}
