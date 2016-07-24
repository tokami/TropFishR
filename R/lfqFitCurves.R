#' @title Fitting of VBGF curves through length frequency data (ELEFAN 1)
#'
#' @description Second step of the Electronic LEngth Frequency ANalysis (ELEFAN), which fits von Bertalanffy
#'    growth curves through restructured length-frequency data.
#'
#' @param lfq a list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column);
#' }
#' @param par a list with following growth parameters:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default: 100),
#'   \item \strong{K} curving coefficient (default: 0.1),
#'   \item \strong{C} amplitude of growth oscillation (default: 0),
#'   \item \strong{WP} winter point (WP = ts + 0.5) (default: 0);
#' }
#'
#' @examples
#' data(trout)
#' res <- lfqRestructure(trout)
#' lfqFitCurves(res)
#'
#' @details This function is used in the analysis of growth parameters with the \code{\link{ELEFAN}} function. It is often
#'    referred to as ELEFAN 1. C expresses the amplitude of growth oscillations and
#'    should be between 0 (no oscillation) and 1 (oscillation with period of zero growth), values above 1 suggest
#'    a period of negative growth, and thus are unlikely. Winter point (WP) designates the period of the year when
#'    growth is slowest, which
#'    in the northern hemisphere is often found around 0.2 (February) and in the southern
#'    hemisphere aorund 0.7 (October) (Pauly and Morgan, 1987).
#'
#' @return A list with following list objects:
#' \itemize{
#'   \item \strong{max_ESP}: maximum explained sum of peaks,
#'   \item \strong{startingSample}: sample number which yield in best ESP value,
#'   \item \strong{startingLength}: length (in cm) which yield in best ESP value;
#' }
#'
#' @importFrom stats na.omit
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a revised and expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. 1981. The relationship between gill surface area and growth performance in fish:
#' a generalization of von Bertalanffy's theory of growth. \emph{Meeresforsch}. 28:205-211
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Pauly, D., 1985. On improving operation and use of ELEFAN programs. Part I: Avoiding
#' "drift" of K towards low values. \emph{ICLARM Conf. Proc.}, 13-14
#'
#' Pauly, D., 1987. A review of the ELEFAN system for analysis of length-frequency data in
#' fish and aquatic invertebrates. \emph{ICLARM Conf. Proc.}, (13):7-34
#'
#' Pauly, D. and G. R. Morgan (Eds.), 1987. Length-based methods in fisheries research.
#' (No. 13). WorldFish
#'
#' Pauly, D. and G. Gaschuetz. 1979. A simple method for fitting oscillating length growth data, with a
#' program for pocket calculators. I.C.E.S. CM 1979/6:24. Demersal Fish Cttee, 26 p.
#'
#' Pauly, D. 1984. Fish population dynamics in tropical waters: a manual for use with programmable
#' calculators (Vol. 8). WorldFish.
#'
#' Quenouille, M. H., 1956. Notes on bias in estimation. \emph{Biometrika}, 43:353-360
#'
#' Somers, I. F., 1988. On a seasonally oscillating growth function. ICLARM Fishbyte 6(1): 8-11.
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2): 407 p.
#'
#' Tukey, J., 1958. Bias and confidence in not quite large samples.
#' \emph{Annals of Mathematical Statistics}, 29: 614
#'
#' Tukey, J., 1986. The future of processes of data analysis. In L. V. Jones (Eds.),
#' The Collected Works of John W. Tukey-philosophy and principles of data analysis:
#' 1965-1986 (Vol. 4, pp. 517-549). Monterey, CA, USA: Wadsworth & Brooks/Cole
#'
#' @export

lfqFitCurves <- function(lfq, par = list(Linf = 100, K = 0.1, C = 0, WP = 0)){

  classes <- lfq$midLengths
  interval <- (classes[2] - classes[1])/2
  n_classes <- length(classes)
  catch <- lfq$catch
  n_samples <- dim(catch)[2]
  peaks_mat <- lfq$peaks_mat
  catch_aAF_F <- lfq$rcounts
  sp.years <- lfq$samplingPeriod
  days.years <- lfq$samplingDays
  delta_list <- lfq$delta_t[,,1]


  # Linf and K for this run
  if(!("Linf" %in% names(par)) | !("K" %in% names(par))) stop("At least 'Linf' and 'K' have to be defined in par!")
  Linfi <- get("Linf",par)
  Ki <- get("K",par)
  if("C" %in% names(par)){
    C <- get("C",par)
  }else C <- 0
  if("WP" %in% names(par)){
    WP <- get("WP",par)
  }else WP <- 0
  if("ts" %in% names(par)) WP <- 0.5 + get("ts",par)


  # maximum length
  Lmax <-  Linfi * 0.98

  # maximum age
  tmax <-  log(1-(Lmax/Linfi)) / - Ki


  # lookup table for soVBGF
  if(C != 0 | WP != 0){
    lookup_age <- seq(0,(tmax+30),0.01)  # initial tmax guess (without soVBGF) plus 10 to be on the save side!
    ## OLD: lookup_length <- Linfi * (1 - exp(-Ki * (lookup_age) - ((C*Ki)/(2*pi)) *((sin(2*pi*(lookup_age-(WP-0.5)))) - (sin(2*pi*(-(WP-0.5)))))))
    lookup_length <- Linfi * (1 - exp(-Ki * (lookup_age-0) +
                                        (((C*Ki)/(2*pi)) * sin(2*pi*(lookup_age-(WP-0.5)))) -
                                        (((C*Ki)/(2*pi)) * sin(2*pi*(0-(WP-0.5))))))

  }
  #plot(lookup_age,lookup_length,type = 'l')

  # in case of soVBGF
  if(C != 0 | WP != 0){
    lookup_ind <- which.min(abs(lookup_length - Lmax))
    tmax <- lookup_age[lookup_ind][1]
  }


  lower_classes <- classes - interval
  upper_classes <- classes + interval


  # starting points meaning combination of sampling time and lengths (in Fisat by 0.5 cm)
  starting_lengths <- seq((classes[1]-interval),(classes[n_classes]+interval),0.5)
  n_startlengths <- length(starting_lengths)


  pot_ages <- log(1-starting_lengths/Linfi) / -Ki
  # delete everything what is exceeding max age (tmax)
  pot_ages[which(pot_ages > tmax)] <- NA
  pot_birthdays <- -pot_ages

  # how many cohorts? and how many older how many younger? depends on starting sample (and class?)
  # the oldest cohort is always the one dying at upper boundary of largest length class (or Linf at sampling time s1)
  # or which has Lmax or tmax at s1
  # youngest cohort is always the one which is born at last sampling time
  # cohort number should always be the same no matter the starting point because the samples stay the same
  ##### OLD:  n_cohorts <- length(-tmax : sp.years)
  # how many cohorts are younger how many older?
  delta_mat <- do.call(rbind,delta_list)
  birthday_actual_cohort <- matrix(rep(delta_mat[1,],each = n_startlengths),ncol=n_samples,nrow=n_startlengths) + pot_birthdays  #### OLD: birthday_actual_cohort <- pot_birthdays[1,]
  birthday_oldest_cohort <- matrix(rep(-tmax + delta_mat[,1],each = n_startlengths),ncol=n_samples,nrow=n_startlengths) # for later if starting samples are added: + delta_list[[first element of each list element]]
  birthday_youngest_cohort <-  matrix(rep(delta_mat[,n_samples],each = n_startlengths),ncol=n_samples,nrow=n_startlengths) # delta_list[[1]][n_samples]  # for later when different starting samples are added: delta_list[[last element of each list element]]

  older_cohorts <- floor(abs(birthday_oldest_cohort - birthday_actual_cohort))   # to account that the variation potentially creates older cohorts
  younger_cohorts <- floor(abs(birthday_actual_cohort - birthday_youngest_cohort))  # to account that the variation potentially creates younger cohorts

  n_cohorts_pre <- older_cohorts + younger_cohorts
  n_cohorts_cor <- max(n_cohorts_pre, na.rm = TRUE)     # for array n_cohorts cant have different length -> no problem to have more cohorts then necessary just more calculation
  older_cohorts <- n_cohorts_cor - younger_cohorts   # add more older cohorts, because there more younger ones any way
  n_cohorts <- n_cohorts_cor + 1 # plus 1 for the actual cohort (which is 0 in -older_cohorts:younger_cohorts)

  # save the dimensions for creation of all arrays later
  dimensions = c(n_samples, n_cohorts, n_startlengths, n_samples)

  # create array with the younger (+) and older (-) cohorts to add later to the ages
  # in two loops because for each starting point (sampling time and class) there is a different amount of
  # younger and older cohorts
  add_cohort_array <- array(0, dim = dimensions)
  for(sampli in 1:n_samples){
    older_cohorts_i <- older_cohorts[,sampli]
    younger_cohorts_i <- younger_cohorts[,sampli]
    for(classi in 1:n_startlengths){
      if(!is.na(older_cohorts_i[classi]) & !is.na(younger_cohorts_i[classi])){
        temp_vec <- rep(-older_cohorts_i[classi]:younger_cohorts_i[classi], each=n_samples)
        temp_array <- array(temp_vec, dim = c(n_samples, n_cohorts))
        add_cohort_array[,,classi,sampli] <- add_cohort_array[,,classi,sampli] + temp_array
      }
    }
  }

  # VBGF
  if(C == 0 | WP == 0) pot_ages_pre <- (log(1-starting_lengths/Linfi) / -Ki)
  # in case of soVBGF
  if(C != 0 | WP != 0){
    lookup_ind_pre <- lapply(starting_lengths, FUN = function(x) which.min(abs(lookup_length - x))[1])
    lookup_ind_pre <- unlist(lookup_ind_pre)
    pot_ages_pre <- lookup_age[lookup_ind_pre]
  }

  pot_ages_array <- array(rep(rep(pot_ages_pre,each=n_samples*n_cohorts),n_startlengths),
                          dim = dimensions)

  pot_ages_array <- pot_ages_array - add_cohort_array
  pot_birthdays <- -pot_ages_array

  # delete all individuals which would be above tmax
  pot_ages_array[which(pot_ages_array > tmax)] <- NA


  delta_array <- array(NA,dim=dimensions)
  for(sampli in 1:n_samples){
    temp_array <- array(rep(delta_list[[sampli]],n_cohorts*n_startlengths),
                        dim = c(n_samples, n_cohorts, n_startlengths))
    delta_array[,,,sampli] <- temp_array
  }

  ages_all_cohorts <- pot_ages_array + delta_array

  # translate into lenth
  # OLD: length_array <- Linfi * (1 - exp(-Ki * ages_all_cohorts) - ((C*Ki)/(2*pi)) *((sin(2*pi*(ages_all_cohorts-(WP-0.5)))) - (sin(2*pi*(-(WP-0.5))))))
  length_array <- Linfi * (1 - exp(-Ki * ages_all_cohorts + (((C*Ki)/(2*pi)) * sin(2*pi*(ages_all_cohorts-(WP-0.5)))) - (((C*Ki)/(2*pi))*sin(2*pi*(-(WP-0.5))))))

  length_array[which(length_array > Lmax)] <- NA


  # LAST STEP! WORKS VERY FAST AND WORKS ALSO IF ALL COHORTS ARE ESTIMATED
  # which classes are hit?

  length_array[length_array < (lower_classes[1]-0.00001)] <- NA     # plus and minus very small number is necessary because length_mat is calculated at 14 digits
  length_array[length_array > (upper_classes[n_classes]+0.00001)] <- NA

  lower_classes_X <- lower_classes-0.00001
  class_vec <- findInterval(length_array, lower_classes_X)  # right open intervals! left open? also
  class_array <- array(class_vec, dim = dimensions)


  # get from class_array to esp_array
  # get indicator_array which contains numbers which refer to specific position within loop_mat
  indicator_array <- array(rep(seq(0,n_classes*(n_samples-1),n_classes),
                               n_cohorts*n_startlengths*n_samples),
                           dim = dimensions)
  indicator_array <- class_array + indicator_array

  cohort_all_sampling_times_array <- apply(indicator_array, MARGIN = c(3,4), FUN = function(x) na.omit(as.vector(x))) #c(na.omit(c(x))))

  # get asp_mat
  loop_mat <- as.matrix(catch_aAF_F)

  # two rounds? or all in one? but save variation, meaning row number!
  #vary_vec[XXXX] # is the best shift within the class of best smapling time

  #cohort_all_sampling_times_array


  get.best.ESP <- function(x){
    x_vec <- x[[1]] ##  unlist(x)      # unlist x
    ## peaks_unique <- (!duplicated(peaks_mat[x_vec]) & peaks_mat[x_vec] > 0)  | peaks_mat[x_vec] <= 0
    pmx <- peaks_mat[x_vec]
    scoremx <- loop_mat[x_vec]

    # deal with positive ones
    uni <- unique(pmx)
    positives <- rep(NA,length(uni))
    for(i in 1:length(uni)){
      loopi <- uni[i]
      indi <- which(pmx == loopi)
      positives[i] <- max(scoremx[indi])
    }

    # sum up all negative ones
    neg_indi <- which(pmx == 0)
    negatives <- scoremx[neg_indi]

    # combine
    ESP <- c(positives,negatives)

    # OLD : wrong because chooses always the first score value of a cohort peak which is hit more often, should choose the highest
    # peaks_unique <- (!duplicated(pmx) & pmx > 0) | pmx <= 0
    # hits_whole_peaks <- x_vec[peaks_unique]
    # ESP <- loop_mat[hits_whole_peaks]                                      # get negative and unique positive scores from score matrix

    # and add up
    if(length(ESP) > 0){
      ESP_sum <- sum(ESP,na.rm=TRUE)                             # take the sum
    }else ESP_sum <- -Inf                    # if startingLength exceeds Linf it is numeric(0) then fit is 0 but should be smaller tthan this
    return(ESP_sum)
  }


  ESP_mat <- apply(cohort_all_sampling_times_array, MARGIN = c(1,2), FUN = get.best.ESP)


  max_ESP <- max(ESP_mat,na.rm=TRUE)
  position_best_ESP <- which(ESP_mat == max_ESP,arr.ind = T)[1,]   # [1,] in case there are more than one maxima
  startingSample <- as.numeric(position_best_ESP[2])
  startingLength <- starting_lengths[as.numeric(position_best_ESP[1])]

  ret <- list(max_ESP = max_ESP,
              startingSample = startingSample,
              startingLength = startingLength)
  return(ret)
}



