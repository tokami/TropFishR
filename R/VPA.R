#' @title Virtual Population Analysis (VPA)
#'
#' @description This function applies the Virtual Population Analysis (VPA) or
#'    Cohort analysis (CA). Methods used to estimate stock biomass and fishing
#'    mortality per age/length group.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \code{midLengths} or \code{age}: midpoints of the length class
#'   (length-frequency data) or ages (age composition data),
#'   \item \code{Linf}: infinite length for investigated species in cm [cm],
#'   \item \code{K}: growth coefficent for investigated species per year [1/year],
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \code{M}: natural mortality [1/year],
#'   \item \code{a}: length-weight relationship coefficent (W = a * L^b),
#'   \item \code{b}: length-weight relationship coefficent (W = a * L^b),
#'   \item \code{catch}: catch as vector for pseudo cohort analysis,
#'      or a matrix with catches of subsequent years to follow a real cohort;
#' }
#' @param terminalF a fishing mortality value which is used as the terminal FM for the
#'   last age/length group.
#' @param analysis_type determines which type of assessment should be done,
#'   options: "VPA" for classical age-based VPA, "CA" for age- or length-based
#'   Cohort analysis
#' @param catch_corFac optional; correction factor for catch, in case provided
#'   catch does spatially or temporarily not reflect catch for fishing ground of
#'   a whole year.
#' @param algorithm an Algorithm to use to solve for fishing mortality. The default
#'   setting \code{"new"} uses \code{\link[stats]{optimize}},
#'   while \code{"old"} uses the algorithm described by Sparre and Venema (1998).
#' @param plot logical; indicating whether a plot should be printed
#'
#' @details The main difference between virtual population analysis (VPA) and cohort
#'    analysis (CA) is the step of calculating the fishing mortality per age class or
#'    length group. While CA works with an approximation by assuming that all fish are
#'    caught during a single day, which makes the calcualtion easier, VPA assumes that
#'    the fish are caught continuously, which has to be solved by the trial and error
#'    method (Sparre and Venema, 1998).
#'    The catch has to be representative for fished species, that means there should not be
#'    other fisheries fishing the same stock. If this is the case \code{catch_corFac} can
#'    be used as a raising factor to account for the proportion of fish caught by other
#'    fisheries.
#'    When the model should follow a real cohort instead of a pseudo cohort, \code{catch}
#'    has to be provided as matrix. The model then starts to follow the first age class
#'    in the first column.
#'    If \code{catch} matrix is shorter than the number of age classes, the age or length
#'    classes without catch information are omitted. It is recommended to only
#'    follow a real cohort if there is enough information for all age classes
#'    (test with: \code{dim(catch)[1] <= dim(catch)[2]}).
#'
#' @keywords function VPA mortality F stock biomass cohort
#'
#' @examples
#' #_______________________________________________
#' # Virtual Popuation Analysis with age-composition data
#' data(whiting)
#' output <- VPA(param = whiting, terminalF = 0.5, analysis_type = "VPA")
#'
#'#_______________________________________________
#' # Pope's Cohort Analysis with age-composition data
#' data(whiting)
#' VPA(whiting, terminalF = 0.5, analysis_type = "CA")
#'
#'#_______________________________________________
#' # Jones's Cohort Analysis with length-composition data
#' data(hake)
#' VPA(hake, terminalF = 0.5, analysis_type = "CA")
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{classes.num}: numeric age classes or length groups (without plus sign),
#'   \item \strong{catch.cohort}: a vector with the catch values which were used for the analysis (exists only if catch was a matrix),
#'   \item \strong{FM_calc}: a vector with the ifshing mortality (M),
#'   \item \strong{Z}: a vector with the total mortality (Z),
#'   \item \strong{survivors}: a vector with the number of fish surviving to the next age class or length group,
#'   \item \strong{annualMeanNr}: ta vector with the mean number of fish per year,
#'   \item \strong{meanBodyWeight}: a vector with the mean body weight in kg,
#'   \item \strong{meanBiomassTon}: a vector with the mean biomass in tons,
#'   \item \strong{YieldTon}: a vector with the yield in tons,
#'   \item \strong{natLoss}: a vector with the number of fish died due to natural mortality,
#'   \item \strong{plot_mat}: matrix with rearranged survivors, nat losses and catches for plotting;
#' }
#'
#'
#' @references
#' Jones, R., 1984. Assessing the effects of changes in exploitation pattern using length
#' composition data (with notes on VPA and cohort analysis). \emph{FAO Fish.Tech.Pap.},
#' (256): 118p.
#'
#' Jones, R., 1990. Length-cohort analysis: the importance of choosing the correct growth
#' parameters. \emph{Journal du Conseil: ICES Journal of Marine Science}, 46(2), 133-139
#'
#' Pope, J.G., 1972. An investigation of the accuracy of virtual population analysis using
#' cohort analysis. \emph{Res.Bull.ICNAF}, (9):65-74
#'
#' Pope, J.G., 1979. A modified cohort analysis in which constant natural mortality is
#' replaced by estimates of predation levels. \emph{ICES C.M.} 1979/H:16:7p. (mimeo)
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' References for weight-length relationship parameters (a & b):
#' Dorel, D., 1986. Poissons del'Atlantique nord-est relations taille-poids.
#' Institut Francais de Recherche pour l'Exploitation de la Mer. Nantes, France. 165 p.
#'
#' @export

VPA <- function(param, terminalF, analysis_type, catch_corFac = NA,
                algorithm="new", plot = FALSE){

  res <- param
  catch <- res$catch

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #           AGE BASED VPA AND COHORT ANALYSIS              #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#

  if("age" %in% names(res) == TRUE){

    classes <- as.character(res$age)

    # Error message if catch and age do not have same length
    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      if(length(classes) != length(catch[,1])) stop("Age/length classes and catch do not have the same length!")
      if(length(classes) != dim(catch)[2]) writeLines("Age/length classes and the real cohort in the catch matrix \ndo not have the same length. The missing age/length \nclasses will be omitted.")
    }else if(class(catch) == 'numeric'){
      if(length(classes) != length(catch)) stop("Age/length classes and catch do not have the same length!")
    }

    a <- res$a
    b <- res$b
    M <- res$M

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    if(class(catch) == 'matrix' | class(catch) == 'data.frame'){
      #find cohort to analyse
      real.cohort <- diag(as.matrix(catch))
       catch.cohort <- c(real.cohort,
                         rep(NA,length(classes.num) - length(real.cohort)))
      if(length(classes.num) != length(real.cohort)){
        catch.cohort <- real.cohort
        classes.num <- classes.num[1:length(catch.cohort)]
      }
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
    lastLengthClass <- max(which(!is.na(catch_cor)),na.rm=TRUE)  ###
    survivors[lastLengthClass] <-
      catch_cor[lastLengthClass] / ((terminalF/(terminalF + M)) * (1 - exp(-(terminalF + M))))


    #   Age-based Cohort Analysis (Pope's cohort analysis)
    if(analysis_type == "CA"){
      # other survivors
      for(x3 in (lastLengthClass-1):1){
        survivors[x3] <- (survivors[x3+1] * exp((M/2)) +
                            catch_cor[x3] ) * exp((M/2))
      }

      #F
      FM_calc <- rep(NA,length(classes.num))
      FM_calc[lastLengthClass] <- terminalF
      for(x5 in 1:(lastLengthClass-1)){
        FM_calc[x5] <- log(survivors[x5]/survivors[x5+1]) - M
      }
    }

    # Traditional VPA
    if(analysis_type == "VPA"){
      #other survivors and fishing mortality
      FM_calc <- rep(NA,length(classes.num))
      FM_calc[lastLengthClass] <- terminalF

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
        FM_calc[num_class] <- sur.F

        #fill survivors
        survivors[num_class] <- survivors[(num_class+1)] *
          exp(sur.F + sur.M)
      }
    }

    # Z
    Z <- rep(NA,length(classes.num))
    for(x6 in 1:(length(Z))){
      Z[x6] <- M  +  FM_calc[x6]
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

    #save all in list
    ret <- c(res,list(
      classes.num = classes.num,
      catch.cohort = catch.cohort,
      FM_calc = FM_calc,
      Z = Z,
      survivors = survivors,
      annualMeanNr = annualMeanNr,
      meanBodyWeight = meanBodyWeight,
      meanBiomassTon = meanBiomassTon,
      yieldTon = yieldTon,
      natLoss = natLoss,
      plot_mat = df.VPAnew))

    class(ret) <- "VPA"

    # plot results
    if(plot==TRUE) try(plot(ret))

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

  if("midLengths" %in% names(res) == TRUE & analysis_type == "VPA") stop("Please choose
                                                         analysis_type = 'CA' for
                                                         length composition data!")

  if((class(catch) == 'numeric' | class(catch) == 'integer') &
     "midLengths" %in% names(res) == TRUE &
     analysis_type == "CA"){

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
    FM_calc <- rep(NA,length(classes.num))
    for(x5 in 1:(length(FM_calc))){
      FM_calc[x5] <- M  *  (F_Z[x5] / (1 - F_Z[x5]))
    }

    # Z
    Z <- rep(NA,length(classes.num))
    for(x6 in 1:(length(Z))){
      Z[x6] <- M  +  FM_calc[x6]
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

    #save all in list
    ret <- c(res,list(
      classes.num = classes.num,
      FM_calc = FM_calc,
      Z = Z,
      survivors = survivors,
      annualMeanNr = annualMeanNr,
      meanBodyWeight = meanBodyWeight,
      meanBiomassTon = meanBiomassTon,
      yieldTon = yieldTon,
      natLoss = natLoss,
      plot_mat = df.VPAnew))

    class(ret) <- "VPA"

    # plot results
    if(plot == TRUE) try(plot(ret))

    return(ret)
  }
}
