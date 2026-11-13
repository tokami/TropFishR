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
#'   \item \code{M}: natural mortality [1/year] (numeric value or vector of identical
#'      length than midLengths),
#'   \item \code{a}: length-weight relationship coefficent (W = a * L^b; in kg/cm^b,
#'      or in g/cm^b with \code{LW_unit = "g"}),
#'   \item \code{b}: length-weight relationship coefficent (W = a * L^b),
#'   \item \code{catch}: total catch of the stock per age class or length
#'      group, as vector for pseudo cohort analysis, or a matrix with catches
#'      of subsequent years to follow a real cohort (or, for length-frequency
#'      data, one column per year selected with \code{catch_columns}).
#'      By default the catch is in numbers (see \code{catch_unit}). For
#'      length-frequency data with the sampled numbers of fish, provide the
#'      total landings with \code{landings} instead (see Details).
#' }
#' @param catch_columns numerical; indicating the column of the catch matrix which should be
#'   used for the analysis.
#' @param catch_unit optional; unit of \code{catch}: \code{NA} (default;
#'    numbers of fish), thousand or million individuals (\code{"'000"},
#'    \code{"'000000"}), or, for the length-based VPA/CA only, the catch in
#'    weight per length group (\code{"kg"} or \code{"tons"}), which is
#'    converted into numbers with the mean body weight of each length group.
#'    Note that this is the weight of the catch of each length group, not
#'    the total landings: for the total landings use \code{landings}.
#' @param catch_corFac optional; additional raising factor for the catch, in
#'   case the catch does spatially or temporarily not reflect the catch of
#'   the stock in a whole year (e.g. landings of only one of several fleets).
#'   With \code{boot}, one value per sample is possible.
#' @param terminalF the fishing mortality rate of the last age/length group.
#' @param terminalE the exploitation rate of the last age/length group.
#' @param analysis_type determines which type of assessment should be done,
#'   options: "VPA" for age or length-based Virtual Population Analysis, "CA" for age- or length-based
#'   Cohort Analysis. Default is "VPA".
#' @param algorithm an Algorithm to use to solve for fishing mortality. The default
#'   setting \code{"new"} uses \code{\link[stats]{optimise}},
#'   while \code{"old"} uses the algorithm described by Sparre and Venema (1998).
#' @param plus_group logical; indicating if the last length group is a plus group (default: TRUE).
#' @param plot logical; indicating whether a plot should be printed
#' @param LW_unit unit of the weight in the length-weight relationship
#'    (\code{a}): \code{"kg"} (default; W in kg) or \code{"g"} (W in g, then
#'    \code{a} is divided by 1000). The mean body weight is in kg and biomass
#'    and yield in tonnes.
#' @param landings optional; total landings (or catch) of the stock in kg,
#'    one value per sample (\code{raiseBy = "sample"}) or year
#'    (\code{raiseBy = "year"}). If provided, \code{param} has to be
#'    length-frequency data with the sampled numbers of fish, which are raised
#'    to the total catch with \code{\link{lfqRaise}} and summed per year
#'    (select the year with \code{catch_columns} if there are several).
#'    With \code{boot}, each resample is raised.
#' @param sampleWeight optional; with \code{landings}, the measured weight of
#'    the samples in kg (one value per sample or year). By default estimated
#'    from the sampled numbers and the length-weight relationship (see
#'    \code{\link{lfqRaise}}).
#' @param raiseBy with \code{landings}, \code{"sample"} (default) to raise
#'    each sample separately or \code{"year"} to raise all samples of a year
#'    with one factor (see \code{\link{lfqRaise}}).
#' @param boot optional; an object of class `lfqBoot` with bootstrapped growth
#'    parameters, natural mortality (see `natMort`) and fishing mortality
#'    (`FM`, from \code{\link{catchCurve}} with `natMort`). If provided, the
#'    length-based VPA (or cohort analysis) is applied to each resampled data
#'    set (re-created with \code{\link{lfqResample}} and the seeds of `boot`)
#'    with the growth parameters and mortality rates of the resample: the
#'    terminal fishing mortality is the `FM` of the resample and a plus group is
#'    created at Linf of the resample. `param` has to be the lfq data used for
#'    the bootstrapped growth analysis and contain the parameters `a` and `b`
#'    of the length-weight relationship.
#' @param natMort with `boot`, the name of the column in `boot$bootRaw` with
#'    natural mortality estimates (e.g. "M_Then").
#' @param yearSel optional; with `boot`, year(s) to be used if the lfq data
#'    covers several years.
#' @param yearCombine logical; with `boot`, should the catches of several years
#'    be combined? Default = FALSE
#' @param binSize optional; with `boot`, new bin size for the resampled data.
#' @param CI with `boot`, percentage for confidence intervals (Default: 95)
#'
#' @details The main difference between virtual population analysis (VPA) and cohort
#'    analysis (CA) is the step of calculating the fishing mortality per age class or
#'    length group. While CA works with an approximation by assuming that all fish are
#'    caught during a single day, which makes the calcualtion easier, VPA assumes that
#'    the fish are caught continuously, which has to be solved by the trial and error
#'    method (Sparre and Venema, 1998).
#'    The VPA/CA requires the total catch of the stock in one year (absolute
#'    numbers), while most other methods (e.g. ELEFAN, \code{\link{catchCurve}})
#'    only use the relative length composition of the samples. For
#'    length-frequency data, keep the sampled numbers of fish in the lfq data
#'    and provide the total landings in weight (\code{landings}): the samples
#'    are then raised to the total catch in numbers by
#'    \code{\link{lfqRaise}}. Alternatively, \code{catch} can contain the
#'    total catch per age class or length group in numbers (or '000 numbers)
#'    or, for the length-based VPA/CA, in weight per length group
#'    (\code{catch_unit}).
#'    The catch has to be representative for fished species, that means there should not be
#'    other fisheries fishing the same stock. If this is the case \code{catch_corFac} can
#'    be used as a raising factor to account for the proportion of fish caught by other
#'    fisheries.
#'    Raising the catch with a single factor changes the numbers and biomass,
#'    but not the fishing mortality at length; different raising factors per
#'    sample change the length composition and, thus, also the fishing
#'    mortality at length.
#'    When the model should follow a real cohort instead of a pseudo cohort, \code{catch}
#'    has to be provided as matrix. The model then starts to follow the first age class
#'    in the first column.
#'    If \code{catch} matrix is shorter than the number of age classes, the age or length
#'    classes without catch information are omitted. It is recommended to only
#'    follow a real cohort if there is enough information for all age classes
#'    (test with: \code{dim(catch)[1] <= dim(catch)[2]}).
#'    If \code{plus_group} is TRUE a different calculation for the survivors of the last length group
#'    is used (for more details please refer to Sparre & Venema (1998)).
#'
#' @keywords function VPA mortality F stock biomass cohort
#'
#' @examples
#' #_______________________________________________
#' # Virtual Popuation Analysis with age-composition data
#' data(whiting)
#' output <- VPA(param = whiting, catch_columns = 1, terminalE = 0.5, analysis_type = "VPA")
#' plot(output)
#'#_______________________________________________
#' # Pope's Cohort Analysis with age-composition data
#' data(whiting)
#' VPA(whiting, terminalE = 0.5, catch_columns = 3, analysis_type = "CA",
#'    plot= TRUE, plus_group = TRUE)
#'
#'#_______________________________________________
#' # Virtual population analysis with length-composition data
#' data(hake)
#' VPA(hake, terminalE = 0.5, analysis_type = "VPA", plot = TRUE,
#'     catch_unit = "'000", plus_group = TRUE)
#'#_______________________________________________
#' # Jones's Cohort Analysis with length-composition data
#' data(hake)
#' VPA(hake, terminalE = 0.5, analysis_type = "CA", plot = TRUE,
#'    catch_unit = "'000", plus_group = TRUE)
#'
#'#_______________________________________________
#' # Length-based VPA with sampled length-frequency data and total landings
#' data(synLFQ7)
#' lfq <- synLFQ7
#' lfq$par <- list(Linf = 123, K = 0.2, t0 = 0, M = 0.2, a = 0.015, b = 3)
#' ## plus group below Linf
#' pg <- lfq$midLengths[max(which(lfq$midLengths < lfq$par$Linf))]
#' lfq <- lfqModify(lfq, plus_group = pg)
#' ## landings in kg at the sampling dates (one value per sample)
#' landings <- seq(20000, 40000, length.out = length(lfq$dates))
#' res <- VPA(lfq, landings = landings, terminalE = 0.5, LW_unit = "g",
#'            plus_group = TRUE)
#' res$raising
#' plot(res)
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{classes.num}: numeric age classes or length groups (without plus sign),
#'   \item \strong{catch.cohort}: a vector with the catch values which were used for
#'   the analysis (exists only if catch was a matrix),
#'   \item \strong{FM_calc}: a vector with the ifshing mortality (M),
#'   \item \strong{Z}: a vector with the total mortality (Z),
#'   \item \strong{survivors}: a vector with the number of fish surviving to the
#'       next age class or length group (same unit than input catch vector),
#'   \item \strong{annualMeanNr}: ta vector with the mean number of fish per year
#'      (same unit than input catch vector),
#'   \item \strong{meanBodyWeight}: a vector with the mean body weight in kg,
#'   \item \strong{meanBiomassTon}: a vector with the mean biomass in tons,
#'   \item \strong{yieldTon}: a vector with the yield in tons,
#'   \item \strong{natLoss}: a vector with the number of fish died due
#'   to natural mortality,
#'   \item \strong{plot_mat}: matrix with rearranged survivors, nat losses
#'   and catches for plotting,
#'   \item \strong{raising}: with \code{landings}, the raising factors (see
#'   \code{\link{lfqRaise}});
#' }
#' If `boot` is provided, the \code{lfqBoot} object with the additional
#' columns `N` (total mean number of fish) and `B` (total mean biomass in
#' tonnes) in `bootRaw`, `maxDen`, `median` and `CI`. The fishing mortality at
#' length of each resample is stored in `misc$FMvecVPA` (used by
#' \code{\link{predict_mod}} with `boot`). Resamples for which the analysis
#' fails (e.g. FM <= 0) get NA.
#'
#' @importFrom graphics plot
#' @importFrom stats optimise
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

VPA <- function(param,
                catch_columns = NA, catch_unit = NA, catch_corFac = NA,
                terminalF = NA, terminalE = NA,
                analysis_type = "VPA", algorithm = "new",
                plus_group = TRUE, plot = FALSE,
                LW_unit = "kg",
                landings = NULL, sampleWeight = NULL, raiseBy = "sample",
                boot = NULL, natMort = NULL, yearSel = NA,
                yearCombine = FALSE, binSize = NA, CI = 95){

    ## VPA with bootstrapped growth parameters and mortality rates
    if(!is.null(boot)){
        return(VPA_boot(param, boot = boot, catch_unit = catch_unit,
                        catch_corFac = catch_corFac, LW_unit = LW_unit,
                        landings = landings, sampleWeight = sampleWeight,
                        raiseBy = raiseBy,
                        analysis_type = analysis_type, algorithm = algorithm,
                        natMort = natMort, yearSel = yearSel,
                        yearCombine = yearCombine, binSize = binSize, CI = CI))
    }

    ## length-weight relationship in kg
    param <- lwUnitA(param, LW_unit)

    ## raise the sampled numbers of fish to the total catch
    raised <- !is.null(landings)
    if(raised){
        if(!("midLengths" %in% names(param)))
            stop("'landings' is only available for the length-based VPA/CA.")
        if(!is.na(catch_unit))
            stop("With 'landings', the catch is raised to numbers of fish: please leave catch_unit = NA.")
        param <- lfqRaise(param, landings = landings, sampleWeight = sampleWeight,
                          by = raiseBy, LW_unit = "kg")
        ## sum the raised catches per year (columns = years)
        if(is.matrix(param$catch) && !is.null(param$dates) &&
           length(param$dates) == ncol(param$catch)){
            cy <- t(rowsum(t(param$catch), format(param$dates, "%Y")))
            param$dates <- as.Date(paste0(colnames(cy), "-01-01"))
            param$catch <- if(ncol(cy) == 1) as.numeric(cy) else unname(cy)
        }
    }
    res <- param

    if(is.na(catch_columns[1])) catch <- res$catch
    if(!is.na(catch_columns[1])){
        catchmat <- res$catch[,(catch_columns)]
        if(length(catch_columns) > 1){
            catch <- rowSums(catchmat, na.rm = TRUE)
        }else catch <- catchmat
    }


    ## PARAMETERS
    if("par" %in% names(res)){
        M <- res$par$M
        Linf <- res$par$Linf
        K <- res$par$K
        t0 <- ifelse("t0" %in% names(res$par),res$par$t0,0)
        C <- ifelse("C" %in% names(res$par), res$par$C, 0)
        ts <- ifelse("ts" %in% names(res$par), res$par$ts, 0)
        a <- res$par$a
        b <- res$par$b
    }else{
        M <- res$M
        Linf <- res$Linf
        K <- res$K
        t0 <- ifelse("t0" %in% names(res),res$t0,0)
        C <- ifelse("C" %in% names(res), res$C, 0)
        ts <- ifelse("ts" %in% names(res), res$ts, 0)
        a <- res$a
        b <- res$b
    }
    ## maybe growth parameters in par but not others...
    if(is.null(a)) a <- res$a
    if(is.null(b)) b <- res$b
    if(is.null(M)) M <- res$M


    ## CHECKS
    if(is.null(a) || is.null(b)) stop("VPA requires information about the length-weight relationship. Please provide 'a' and 'b' estimates in res.")
    if(is.null(M)) stop("Please provide a natural mortality estimate 'M' in res.")


  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #           AGE BASED VPA AND COHORT ANALYSIS              #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#

  if("age" %in% names(res) == TRUE){

    classes <- as.character(res$age)

    # Error message if catch and age do not have same length
    if(inherits(catch,'matrix') || inherits(catch,'data.frame')){
      #if(length(classes) != length(catch[,1])) stop("Age/length classes and catch do not have the same length!")
      if(length(classes) != length(diag(as.matrix(catch)))) warning("Age/length classes and the real cohort in the catch matrix do not have the same length. The missing age/length classes will be omitted.")
      }else if(inherits(catch,'numeric')){
      if(length(classes) != length(catch)) stop("Age/length classes and catch do not have the same length!")
    }

    if(length(M) == length(classes)){
      M_vec <- M
    }else if(length(M) > 1){
      writeLines(noquote("The number of natural mortality values does not correspond to the number of length classes. \nOnly the first value will be used."))
      M_vec <- rep(M[1],length(classes))
    }else{
      M_vec <- rep(M, length(classes))
    }
    if(!is.na(terminalF)){
      terminalE <- terminalF / (terminalF + M_vec[length(M_vec)])
    }else if(!is.na(terminalE)){
      terminalF <- M_vec[length(M_vec)] * terminalE / (1 - terminalE)
    }else{
      stop("Please provide either the terminal exploitation rate (terminalE) or the terminal fishing mortality rate (terminalF)!")
    }
    terminalZ <- terminalF + M_vec[length(M_vec)]

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

      if(inherits(catch,"matrix") || inherits(catch,"data.frame")){
      writeLines(noquote("A catch matrix was provided: The VPA/CA will follow the 'real' cohort, assuming yearly intervals in the catch matrix. If you want to perform the VPA/CA with a pseudo cohort please use the argument 'catch_columns' specifying which column(s) of the catch matrix to use."))
      #find cohort to analyse
      real.cohort <- diag(as.matrix(catch))
       catch.cohort <- c(real.cohort,
                         rep(NA,length(classes.num) - length(real.cohort)))
      if(length(classes.num) != length(real.cohort)){
        catch.cohort <- real.cohort
        classes.num <- classes.num[1:length(catch.cohort)]
      }
       if(length(M_vec) != length(classes.num)){
         M_vec <- M_vec[1:length(classes.num)]
       }
    }
    if(inherits(catch,'numeric')){
      catch.cohort <- catch
    }

    #Correct catch  if not representative for one year
    if(!is.na(catch_corFac)) catch_cor <- catch.cohort * catch_corFac
    if(is.na(catch_corFac)) catch_cor <- catch.cohort

    # translate catch into individuals
    if(catch_unit %in% c("tons", "t", "T", "Tons", "tonnes", "Tonnes")){
      stop("The age-based VPA/CA is not applicable to catch in weight. Please set the argument 'catch_unit' either to NA, '000, or '000000 for catch in numbers.")
    }else if(catch_unit %in% c("kg", "Kg", "KG", "kilo", "KILO", "kilogramm", "Kilogramm")){
      stop("The age-based VPA/CA is not applicable to catch in weight. Please set the argument 'catch_unit' either to NA, '000, or '000000 for catch in numbers.")
    }else if(catch_unit %in% c("'000","1000","1e3")){
      catch_numbers <- catch_cor * 1000
    }else if(catch_unit %in% c("'000000","1000000","1e6","'000.000")){
      catch_numbers <- catch_cor * 1000000
    }else if(!is.na(catch_unit)){
      stop(paste0(catch_unit, " not known. Please use either NA, '000, or '000000 for catch in numbers."))
    }else{
      if(!raised) warning("You did not specify catch_unit. The Method assumes that catch is provided in numbers!")
      catch_numbers <- catch_cor
    }

    #Survivors    #N(L1)=(N(L2)*H(L1,L2)+C(L1,L2)) *H(L1,L2)
    survivors <- rep(NA,length(classes.num))

    # survivors last size class
    lastLengthClass <- max(which(!is.na(catch_numbers)),na.rm=TRUE)  ###
    if(!plus_group) survivors[lastLengthClass] <-
      catch_numbers[lastLengthClass] / (terminalE * (1 - exp(-terminalZ)))
    if(plus_group) survivors[lastLengthClass] <- catch_numbers[lastLengthClass] / terminalE

    #   Age-based Cohort Analysis (Pope's cohort analysis)
    if(analysis_type == "CA"){
      # other survivors
      for(x3 in (lastLengthClass-1):1){
        survivors[x3] <- (survivors[x3+1] * exp((M_vec[x3]/2)) +
                            catch_numbers[x3] ) * exp((M_vec[x3]/2))
      }

      #F
      FM_calc <- rep(NA,length(classes.num))
      FM_calc[lastLengthClass] <- terminalF
      for(x5 in 1:(lastLengthClass-1)){
        FM_calc[x5] <- log(survivors[x5]/survivors[x5+1]) - M_vec[x5]
      }
    }

    # Traditional VPA
    if(analysis_type == "VPA"){
      #other survivors and fishing mortality
      FM_calc <- rep(NA,length(classes.num))
      FM_calc[lastLengthClass] <- terminalF

      for(num_class in (lastLengthClass-1):1){

        sur.C <- catch_numbers[num_class]
        sur.Ntplus1 <- survivors[(num_class+1)]
        sur.M <- M_vec[num_class]
        LHS <-  sur.C / sur.Ntplus1
        sur.F <- 0

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
          tmp <- optimise(Fcalc, interval=c(0,100))
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
    Z <- M_vec + FM_calc

    #Annual mean Nr
    deads <- abs(diff(survivors))
    annualMeanNr <- deads / Z[-length(survivors)]
    # annualMeanNr[length(survivors)] <- NA
    annualMeanNr[length(survivors)] <- survivors[length(survivors)] / Z[length(survivors)]

    #FOR PLOT
    #Survivors rearranged
    survivors_rea <- rep(NA,length(classes.num))
    for(x8 in 1:(length(survivors_rea)-1)){
      survivors_rea[x8] <- survivors[x8+1]
    }
    survivors_rea[length(survivors_rea)] <- 0


    #Calculate natural losses
    natLoss <- survivors - survivors_rea - catch_numbers


    #put together in dataframe
    df.VPAnew <- data.frame(survivors = survivors_rea,
                            nat.losses = natLoss,
                            catch = catch_numbers,
                            FM_calc = FM_calc)

    #transpose matrix for barplot function
    df.VPAnew <- t(as.matrix(df.VPAnew))
    colnames(df.VPAnew) <- classes.num

    #save all in list
    ret <- c(res,list(
      classes.num = classes.num,
      catch.cohort = catch.cohort,
      FM_calc = FM_calc,
      Z = Z,
      survivors_L1 = survivors,
      survivors_L2 = survivors_rea,
      annualMeanNr = annualMeanNr,
      natLoss = natLoss,
      plot_mat = df.VPAnew))

    class(ret) <- "VPA"

    # plot results
    if(plot==TRUE) try(plot(ret))

    return(ret)
  }

  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
  #          LENGTH BASED VPA AND COHORT ANALYSIS            #
  #HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH#
    if("midLengths" %in% names(res) == TRUE &&
      (inherits(catch,"matrix") || inherits(catch,"data.frame"))){
    if(is.na(catch_columns[1])) stop("The length-based Cohort analysis is not applicable to length frequency data. Please provide catch as vector or use the argument 'catch_columns'.")
    if(!is.na(catch_columns[1])){
      catch <- rowSums(catch)
    }
  }

  if((inherits(catch,'numeric') || inherits(catch,'integer')) &&
     "midLengths" %in% names(res) == TRUE){

    classes <- as.character(res$midLengths)

      # Error message if catch and age do not have same length
      if(inherits(catch,"matrix") || inherits(catch,"data.frame")){
      if(length(classes) != length(catch[,1])) stop("Midlengths and catch do not have the same length!")
    }else if(inherits(catch,'numeric')){
      if(length(classes) != length(catch)) stop("Midlengths and catch do not have the same length!")
    }

    if(length(M) == length(classes)){
      M_vec <- M
    }else if(length(M) > 1){
      writeLines(noquote("The number of natural mortality values does not correspond to the number of length classes. \nOnly the first value will be used."))
      M_vec <- rep(M[1],length(classes))
    }else{
      M_vec <- rep(M, length(classes))
    }
    if(!is.na(terminalF)){
      terminalE <- terminalF / (terminalF + M_vec[length(M_vec)])
    }else if(!is.na(terminalE)){
      terminalF <- M_vec[length(M_vec)] * terminalE / (1 - terminalE)
    }else{
      stop("Please provide either the terminal exploitation rate (terminalE) or the terminal fishing mortality rate (terminalF)!")
    }
    terminalZ <- terminalF + M_vec[length(M_vec)]

    # correct catch with raising factor
    if(!is.na(catch_corFac)) catch_cor <- catch * catch_corFac
    if(is.na(catch_corFac)) catch_cor <- catch

    # create column without plus group (sign) if present
    classes.num <- do.call(rbind, strsplit(classes, split="\\+"))
    classes.num <- as.numeric(classes.num[,1])

    #calculate size class interval
    interval <- classes.num[2] - classes.num[1]

    # lower and upper length vectors
    lowerLength <- classes.num - (interval / 2)
    upperLength <- classes.num + (interval / 2)
    if(plus_group) upperLength[length(upperLength)] <- Linf

    #Mean body weight
    # FAO manual:
    meanBodyWeight <- a * ((lowerLength + upperLength)/2)^b    # a * classes.num ^ b
    # same as what provided in FAO manual: a * ((lowerLength + upperLength)/2)^b
    #meanBodyWeight <- meanBodyWeight / 1000  # in kg
    #according to Beyer (1987) (FISAT II)
    # meanBodyWeight <- (1/(upperLength - lowerLength)) * (a / (b + 1)) * (upperLength^(b+1) - lowerLength^(b+1))

    # translate catch in tons into numbers
    if(catch_unit %in% c("tons", "t", "T", "Tons", "tonnes", "Tonnes")){
      catch_numbers <- (catch_cor * 1000) / meanBodyWeight
    }else if(catch_unit %in% c("kg", "Kg", "KG", "kilo", "KILO", "kilogramm", "Kilogramm")){
      catch_numbers <- catch_cor / meanBodyWeight
    }else if(catch_unit %in% c("'000","1000","1e3")){
      catch_numbers <- catch_cor * 1000
    }else if(catch_unit %in% c("'000000","1000000","1e6","'000.000")){
      catch_numbers <- catch_cor * 1000000
    }else if(!is.na(catch_unit)){
      stop(paste0(catch_unit, " not known. Please use either 'tons' or 'kg' for catch in weight or NA, '000, or '000000 for catch in numbers."))
    }else{
      if(!raised) warning("You did not specify catch_unit. The Method assumes that catch is provided in numbers!")
      catch_numbers <- catch_cor
    }

    # t of lower length classes
    t_L1 <- t0 - (1/K) * log(1 - (lowerLength / Linf))

    # t of lower upper classes
    t_L2 <- t0 - (1/K) * log(1 - (upperLength / Linf))
    if(upperLength[length(upperLength)] > Linf){
      writeLines(noquote("Upper limit of last length class is larger than Linf, \nconsider creating lower plus group or set the argument plus_group = TRUE."))
    }

    # delta t
    # dt <- t_L2 - t_L1
    dt <- (1/K) * log((Linf - lowerLength)/(Linf - upperLength))

    #Survivors
    survivors <- rep(NA, length(classes.num))

    # survivors last size class
    lastLengthClass <- max(which(!is.na(catch_numbers)),na.rm=TRUE)
    if(plus_group) survivors[length(survivors)] <- catch_numbers[length(survivors)] / terminalE
    if(!plus_group){
      survivors[length(survivors)] <- catch_numbers[length(survivors)] /
        (terminalE * (1 - exp(-terminalZ * dt[length(survivors)])))
    }

    ###  Jones' Length-based Cohort Analysis
    if(analysis_type == "CA"){
      # H (L1,L2)   #H(L1,L2)=((Linf-L1)/Linf-L2)^(M/2K)
      H <- ((Linf - lowerLength)/(Linf - upperLength))^(M_vec/(2*K))

      # other survivors
      for(x3 in (length(survivors)-1):1){
        survivors[x3] <- (survivors[x3+1] * H[x3] + catch_numbers[x3]) * H[x3]
      }

      # F/Z  #F(L1,L2)/Z(L1,L2)=C(L1,L2)/(N(L1)-N(L2))
      deads <- abs(diff(survivors))
      F_Z <- catch_numbers[-length(survivors)] / deads
      F_Z[length(survivors)] <- terminalE

      #F  # F = M * (F_Z / 1-F_Z)
      FM_calc <- M_vec * F_Z / (1 - F_Z)
    }

    ###  Length-based VPA
    if(analysis_type == "VPA"){
      #other survivors and fishing mortality
      FM_calc <- rep(NA,length(classes.num))
      FM_calc[lastLengthClass] <- terminalF

      for(num_class in (lastLengthClass-1):1){

        sur.C <- catch_numbers[num_class]
        sur.Ntplus1 <- survivors[(num_class+1)]
        sur.M <- M_vec[num_class]
        sur.dt <- dt[num_class]
        LHS <-  sur.C / sur.Ntplus1
        sur.F <- 0

        if(algorithm == "old"){
          LHS <-  sur.C / sur.Ntplus1
          sur.F <- 0
          seqi <- c(1e-1,1e-2,1e-3,1e-4,1e-5,1e-6,1e-7)
          #trail and error
          for(y in seqi){
            stepi <- y
            for(x in seq(sur.F,10,stepi)){
              sur.F <- x
              RHS <- (sur.F/(sur.F + sur.M)) * (exp(sur.F+sur.M) * sur.dt - 1)
              if(LHS-RHS < 0) break
            }
            sur.F = x-stepi
          }
        }

        if(algorithm == "new"){
          Fcalc <- function(sur.F=sur.M){
            ((sur.F/(sur.F+sur.M)) * (exp((sur.F+sur.M) * sur.dt) - 1) - (sur.C / sur.Ntplus1))^2
          }
          tmp <- optimise(f = Fcalc, interval=c(0,100))
          sur.F <- tmp$min
        }

        #fill F
        FM_calc[num_class] <- sur.F

        #fill survivors
        survivors[num_class] <- survivors[(num_class+1)] *
          exp((sur.F + sur.M) * sur.dt)
      }
    }

    # Z
    Z <- M_vec + FM_calc

    #Annual mean Nr
    deads <- abs(diff(survivors))
    annualMeanNr <- deads / Z[-length(survivors)]
    # annualMeanNr[length(survivors)] <- NA
    annualMeanNr[length(survivors)] <- survivors[length(survivors)] / Z[length(survivors)]

    #Mean biomass
    meanBiomass <- annualMeanNr * meanBodyWeight
    meanBiomassTon <- meanBiomass / 1000

    #Yield
    yield <- catch_numbers * meanBodyWeight
    yieldTon <- yield / 1000

    #FOR PLOT
    #Survivors rearranged
    survivors_rea <- rep(NA,length(classes.num))
    for(x8 in 1:(length(survivors_rea)-1)){
      survivors_rea[x8] <- survivors[x8+1]
    }
    survivors_rea[length(survivors_rea)] <- 0

    #Calculate natural losses
    natLoss <- survivors - survivors_rea - catch_numbers

    #put together in dataframe
    df.VPAnew <- data.frame(survivors = survivors_rea,
                            nat.losses = natLoss,
                            catch = catch_numbers,
                            FM_calc = FM_calc,
                            meanBodyWeight = meanBodyWeight,
                            meanBiomassTon = meanBiomassTon)

    #transpose matrix for barplot function
    df.VPAnew <- t(as.matrix(df.VPAnew))
    colnames(df.VPAnew) <- classes.num

    #save all in list
    ret <- c(res,list(
      classes.num = classes.num,
      FM_calc = FM_calc,
      Z = Z,
      meanBodyWeight = meanBodyWeight,
      survivors_L1 = survivors,
      survivors_L2 = survivors_rea,
      catch_numbers = catch_numbers,
      annualMeanNr = annualMeanNr,

      meanBiomassTon = meanBiomassTon,

      yieldTon = yieldTon,
      natLoss = natLoss,
      plot_mat = df.VPAnew))

    class(ret) <- "VPA"

    # plot results
    if(plot == TRUE & all(!is.na(survivors)) & all(!is.na(natLoss))) try(plot(ret))

    return(ret)

  }# stop("Please choose analysis_type = 'CA' for length composition data!")
}


## Convert the parameter a of the length-weight relationship (in param or
## param$par) from g to kg, if LW_unit = "g"
lwUnitA <- function(param, LW_unit = "kg"){
    LW_unit <- match.arg(LW_unit, c("kg", "g"))
    if(LW_unit == "g"){
        if(!is.null(param$a)) param$a <- param$a / 1000
        if(is.list(param$par) && !is.null(param$par$a)) param$par$a <- param$par$a / 1000
    }
    param
}
