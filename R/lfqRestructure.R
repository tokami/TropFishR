#' @title Restructuring of length frequency data (ELEFAN 0)
#'
#' @description First step of the Electronic LEngth Frequency ANalysis (ELEFAN), which is restructuring
#'    length-frequency data (lfq). This is done according to a certain protocol, described by many authors (see
#'    References for more information).
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column);
#' }
#' @param MA number indicating over how many length classes the moving average should be performed (default: 5)
#' @param addl.sqrt additional squareroot transformation of positive values according to Brey et al. (1988)
#'    (default: FALSE).
#'    Particularly useful if many observations have a low frequency (<10).
#'
#' @examples
#' data(trout)
#' lfqRestructure(param = trout)
#'
#' @details This function is used in the analysis of growth parameters with the \link{ELEFAN} function. It is often
#'    referred to as ELEFAN 0.
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

lfqRestructure <- function(param, MA = 5, addl.sqrt = FALSE){
  classes <- param$midLengths
  catch <- param$catch

  interval <- (classes[2]-classes[1])/2
  n_samples <- dim(catch)[2]
  n_classes <- length(classes)

  if("dates" %in% names(param)){
    dates.all <- as.Date(param$dates)
    # get length of smapling period
    # continuous time in years
    julian_days <- as.numeric(format(dates.all, format="%Y")) + as.numeric(format(dates.all, format="%j"))/366
    days.years <- julian_days - julian_days[1]  # OLD: #days <- as.numeric(dates.all - as.Date((dates.all[1]))) #days.years <- days/365

    # sampling period # OLD: sample.period.days <- days[length(days)] - days[1]  sp.years <- sample.period.days/365
    sp.years <- days.years[length(days.years)] - days.years[1]

    time_diff_year <- as.numeric(diff(dates.all)/365)
    cum_diff_year <- cumsum(time_diff_year)

    # from dates of sampling times to relative delta ts
    delta_t1 <- c(0,cum_diff_year)
    delta_list <- vector("list",n_samples)
    delta_list[[1]] <- delta_t1
    for(i in 1:(n_samples-1)){
      delta_list[[i+1]] <- delta_t1 - cum_diff_year[i]
    }
    param$samplingPeriod <- sp.years
    param$samplingDays <- days.years
    param$delta_t <- array(delta_list,dim = c(1,length(delta_list),1))
  }

  # A
  # moving average
  # Computation of the moving average (ma) over 5 class intervals; the computation also considers 2 class intervals above and below the highest and lowest recorded length respectively
  pm <- (MA-1)/2 # plus minus
  catch_MA <- matrix(ncol=n_samples, nrow=n_classes, NA)
  add_zero_mat <- matrix(ncol=n_samples, nrow=pm, 0)
  colnames(add_zero_mat) <- colnames(catch)
  enlarge_mat <- rbind(add_zero_mat,catch,add_zero_mat)

  moving_average <- function(x){
    res <- rep(NA,n_classes)
    for(j in (pm+1):(n_classes+pm)){
      res[j-pm] <- sum(x[(j - pm):(j + pm)], na.rm = TRUE) / MA
    }
    return(res)
  }

  catch_MA <- apply(enlarge_mat, MARGIN = 2, FUN = moving_average)

  # B
  # adjusted frequency
  # Computation of quotients Ni/mai. Also, the mean quotient (m') is computed using
  catch_AF <- catch/catch_MA
  #catch_AF[is.na(catch_AF)] <- 0
  catch_AF[catch_AF == Inf] <-  NaN
  mean_quotient <- apply(catch_AF, MARGIN = c(2), FUN = function(x) mean(x,na.rm=TRUE))#sum(x, na.rm = TRUE)/n_classes ) # ISSUE: maybe: mean_quotient per sample or mean quotient per all samples

  # C
  # relative adjusted frequency
  # Division of (Ni/mai) by m' and subtraction of 1 (one) from quotient
  mean_quotient_mat <- matrix(rep(mean_quotient,each = n_classes),ncol = n_samples, nrow = n_classes)  # necessary because otherwise vector is repeated for length of matrix but not each!
  catch_rAF <- (catch_AF/mean_quotient_mat) - 1

  # Adjustements
  # D
  # Identifying "isolated peaks", i.e. counting in column C the occurrences of zero-frequencies in the two class intervals above and below a given positive value and defining nz = (the number of neighboring classes with zero-frequencies).
  #account for inflation of adjacent zero frequencies = nz
  #low values
  # -1 to 0
  isolated_peaks_mat <- matrix(ncol=n_samples, nrow=n_classes,NA)
  zero_mat <- catch == 0
  add_TF_mat <- matrix(ncol=n_samples, nrow=2,TRUE)
  colnames(add_TF_mat) <- colnames(zero_mat)
  enlarge_zero_mat <- rbind(add_TF_mat, zero_mat, add_TF_mat)

  zero_neighbours <- function(x){
    res <- rep(NA,n_classes)
    for(j in 3:(n_classes+2)){
      res[j-2] <- sum(x[(j-2):(j+2)][-3], na.rm = TRUE)  # minus 3 so that TRUE or FALSE of observed length class is not respected
    }
    return(res)
  }

  isolated_peaks_mat <- apply(enlarge_zero_mat, MARGIN = 2, FUN = zero_neighbours)
  isolated_peaks_mat <- ifelse(catch_rAF>=0, isolated_peaks_mat, isolated_peaks_mat*(-1))  # ISSUE: multiply by minus 1 if in last step smaller zero according to Aaron Greenbarg

  not_effective_mat <- catch == 0 | catch_rAF < 0

  # E
  # De-emphasizing positive values in column C by multiplying them by 0.5nz (rather than 0.2nz in Brey and Pauly 1986).
  catch_aAF <- catch_rAF * (0.5^isolated_peaks_mat)
  catch_aAF[not_effective_mat] <- catch_rAF[not_effective_mat]

  # additional step mentioned in Brey et al. 1988 useful if the data include a lot of low values (<10)
  if(addl.sqrt){
    posFs <- which(catch_aAF > 0)
    if(length(posFs)>0) {catch_aAF[posFs] <- catch_aAF[posFs] / sqrt(1+2/catch[posFs])}
  }

  # F
  # All values in column E that are equal to -1.000 are set to zero. All values that are smaller than zero but not equal to -1.000 are multiplied by SPV/-SNV, where SPV = sum of positive values in column E, and SNV = sum of negative values in column E.
  catch_aAF_F <- catch_aAF
  catch_aAF_F[catch_aAF_F == -1 | is.na(catch_aAF_F)] <- 0

  SPV <- apply(catch_aAF_F, MARGIN = c(2), FUN = function(x) sum(x[x > 0] ,na.rm = TRUE))
  SNV <- apply(catch_aAF_F, MARGIN = c(2), FUN = function(x) abs(sum(x[x < 0] ,na.rm = TRUE)))

  for(i in 1:n_samples){
    TorF2 <- catch_aAF_F[,i] < 0 & catch_aAF_F[,i] != -1 & !is.na(catch_aAF_F[,i]) & !is.nan(catch_aAF_F[,i])
    catch_aAF_F[TorF2,i] <- catch_aAF_F[TorF2,i] * (SPV[i]/SNV[i])
  }


  # another addition, not sure where it comes from....
  #low values
  # TorF <- catch <= 10 & catch > 0 & catch_rAF > 0                    ### ISSUE: not sure if this should be incorporated
  # catch_aAF[TorF] <- catch_aAF[TorF] / (sqrt(1+(2/catch[TorF])))


  # ISSUE: maybe the last two length class (Lmax and Lmax-1) have to be dependent on Linf. because length classes above Linf are not used in the analysis, because log is not defined for them -> NA
  # in addition: decrease influence of largest length groups (p.62 top in Fisat manual '97)
  #ultimate negative to zero
  catch_aAF_F[n_classes,][catch_aAF_F[n_classes,] < 0] <- 0   # positive not modified but "counted only ones" does this apply and is accounted for by flagging out?

  #penultimate negative dividing by 2
  withoutNA_1 <- !is.na(catch_aAF_F[(n_classes-1),])
  withoutNA <- !is.na(catch_aAF_F[(n_classes),])
  # catch_aAF_F[(n_classes-1),][withoutNA & (catch_aAF_F[(n_classes-1),] < 0)] <-
  #   catch_aAF_F[(n_classes-1),][withoutNA & (catch_aAF_F[(n_classes-1),] < 0)] * 0.5

  negPen <- (catch_aAF_F[(n_classes-1),] < 0)
  # posPen <- (catch_aAF_F[(n_classes-1),] >= 0)

  catch_aAF_F[(n_classes-1),][withoutNA_1 & negPen] <-
    catch_aAF_F[(n_classes-1),][withoutNA_1 & negPen] * 0.5   # ISSUE: according to Greenbarg # either copy and divide by 2 the unltimate or penultimate length class

  # catch_aAF_F[(n_classes-1),][withoutNA & posPen] <-
  #   catch_aAF_F[(n_classes),][withoutNA & posPen]  # ISSUE: according to Greenbarg


  param$rcounts <- catch_aAF_F


  # G
  # Computing the available sum of peaks (ASP) by identifying the highest points in a run of positive point values, bordered on either side by zero or by negative values (troughs). Such run may consist of one length class only.
  peaks <- catch_aAF_F > 0
  ranks_mat <- apply(catch_aAF_F,MARGIN = c(2), FUN = function(x) rank(x, ties.method = "first"))
  ranks_mat <- ranks_mat * peaks

  # create peak matrix
  prep_mat <- catch_aAF_F
  prep_mat <- ifelse(prep_mat > 0,1,0)
  peaks_mat <- matrix(NA,ncol = n_samples, nrow = n_classes)
  for(i in 1:n_samples){
    vec_peaki <- prep_mat[,i]
    rle_val <- rle(vec_peaki)$values
    rle_val[which(rle_val == 1)] <- 1:length(rle_val[which(rle_val == 1)])
    peaks_mat[,i] <- rep(rle_val,rle(vec_peaki)$lengths)
    peaks_mat[(peaks_mat[,i]!= 0 & !is.na(peaks_mat[,i])),i] <-
      peaks_mat[(peaks_mat[,i]!= 0 & !is.na(peaks_mat[,i])),i] + 20*i  # to make avery peak unqiue for flagging out later
  }
  param$peaks_mat <- peaks_mat

  # OLD: but this count the max of all peaks! and this is incorrect as I think, better to count ASP which could be hit by one cohort
  ASP_vec <- rep(NA,n_samples)
  for(j in 1:n_samples){
    peaki <- peaks_mat[,j]
    xxx <- catch_aAF_F[,j]
    uni_peaki <- unique(peaki)
    uni_peaki_without0 <- uni_peaki[uni_peaki>0 & !is.na(uni_peaki)]

    max_loopi <- rep(NA,length(uni_peaki_without0))
    for(i in 1:length(uni_peaki_without0)){
      loopi <- xxx[which(peaki == uni_peaki_without0[i])]
      max_loopi[i] <- max(loopi, na.rm = TRUE)
    }
    ASP_vec[j] <- sum(max_loopi,na.rm=TRUE)
  }

  ASP <- sum(ASP_vec, na.rm = TRUE)
  param$ASP <- ASP

  class(param) <- "lfq"
  return(param)
}




