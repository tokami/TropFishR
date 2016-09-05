#' @title Recruitment patterns
#'
#' @description  This function estimates recrutiment patterns from length-frequency data.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \code{midLengths}: midpoints of the length classes (length-frequency
#'   data),
#'   \item \code{Linf}: infinite length for investigated species in cm [cm],
#'   \item \code{K}: growth coefficent for investigated species per year [1/year],
#'   \item \code{t0}: theoretical time zero, when growth curve crosses length equalling zero,
#'   \item \code{D}: optional; for generalised vBGF
#'   \item \code{C}: optional; intensity of the (sinusoid) growth oscillations of
#'   the model,
#'   \item \code{ts}: optional; onset of the positive phase of the growth oscillation (fraction of year; i.e. Jan 1 = 0, July 1 = 0.5, etc.),
#'   \item \code{catch}: catches, vector or matrix with catches of subsequent sampling
#'   times
#' }
#' @param tsample sampling time corresponding to time when catch was sampled as
#'    fraction of year
#'    (e.g. 0.25 for 1st of April). If catch was sampled several times (catch as matrix)
#'    vector has to be provided with sampling times in equal order.
#' @param catch_column numeric; if catch in param is a matrix, this number indicates
#'    the column of the catch matrix which should be used for the analysis.
#' @param plot logical; indicating whether a plot should be printed
#'
#' @examples
#' # one sample
#' dat <- list(midLengths = seq(2,98,4),
#'                catch = c(0.6,17.6,93,83.2,12.6,0.3,0,0,0,1,17.1,51.4,
#'                26.1,2.2,0.2,4.5,21.6,17.6,3.7,8.7,10.6,6.3,5.6,2.9,0.8),
#'                Linf = 100,
#'                K = 0.5,
#'                t0 = 0)
#' recruitment(param = dat, tsample = 0.25)
#'
#'
#' # several samples
#' data(synLFQ4)
#'
#' # add growth parameters
#' synLFQ4$Linf <- 80
#' synLFQ4$K <- 0.5
#' synLFQ4$t0 <- 0.25
#'
#' # retrieve sampling times from catch matrix
#' s_dates <- as.POSIXlt(synLFQ4$dates, format="%d.%m.%Y")
#'
#' recruitment(param = synLFQ4, tsample = s_dates$yday/365, plot = TRUE)
#'
#'
#' plot(synLFQ4, Fname = "catch",
#'    par = list(Linf = 80, K = 0.5, t_anchor = 0.25, C = 0.75, ts = 0),
#'    ylim = c(0,80))
#'
#' @details
#' This function calculates recruitment patterns of a
#' stock by backward projection onto the length axis of a set of length frequency data
#' using the (special, generalised or seasonalised)
#' von Bertallanfy growth curve (vBGF, Pauly 1982). The method assumes that (i) all fish in a
#' data set grow as described by a single set of growth parameters and (ii) one month out
#' of twelve always has zero recruitment. The second assumption is probably not met, since
#' temperate species may contain more than one month with zero recruitment, while tropical species
#' may have more constant recruitment without months of no recruitment.
#' If t0 is not provided, a relative recruitment pattern will be estimated without
#' specific month values returned in the results. However, an estimate of t0
#' can be obtained by the time
#' lag between peak spawning and recruitment. Several length-frequency data sets
#' can be used to estimate the recrutiment pattern by providing catch as a matrix and
#' setting catch_column to NA (default). Then the fraction per time is calculated for
#' each size class in each sample and then pooled together. For the generalised vBGF, D is
#' required, for the seasonalised vBGF C, ts and D.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{ti}: actual age,
#'   \item \strong{tS_frac}: age at which the length was zero expressed as fraction
#'   of the year,
#'   \item \strong{cor_months}: corresponding months,
#'   \item \strong{months}: numeric months (relative if no t0 is not given),
#'   \item \strong{months_abb}: months (only presented if t0 is given),
#'   \item \strong{all_recruits}: number of recruits per month as matrix if several
#'   length-frequency data sets are provided,
#'   \item \strong{mean_recruits}: (mean) number of recruits per month,
#'   \item \strong{per_recruits}: precentage number of recruits per month.
#' }
#'
#' @importFrom graphics plot
#' @importFrom stats aggregate
#'
#' @references
#' Brey, T., Soriano, M., Pauly, D., 1988. Electronic length frequency analysis. A
#' revised and expanded user's guide to ELEFAN 0, 1 and 2. (Second edition).
#' Berichte aus dem Institut f??r Meereskunde Kiel, No 177, 31p.
#'
#' Moreau, J., & Cuende, F. X., 1991. On improving the resolution of the recruitment
#' patterns of fishes. \emph{Fishbyte}, 9(1), 45-46.
#'
#' Pauly, D., 1982. Studying single-species dynamics in a tropical multispecies context.
#' In Theory and management of tropical fisheries. \emph{ICLARM Conference Proceedings}
#' (Vol. 9, No. 360, pp. 33-70).
#'
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

recruitment <- function(param, tsample, catch_column = NA, plot = FALSE){

  res <- param
  Linf <- res$Linf
  K <- res$K
  t0 <- ifelse("t0" %in% names(res),res$t0, 0)
  D <- ifelse("D" %in% names(res),res$D, 1)
  C <- ifelse("C" %in% names(res),res$C, 0)
  ts <- ifelse("ts" %in% names(res),res$ts, 0)
  classes <- as.character(res$midLengths)
  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  if(is.na(catch_column)) catch <- res$catch
  if(!is.na(catch_column)) catch <- res$catch[,catch_column]

  if(class(catch) == "matrix" | class(catch) == "data.frame"){
    len_catch <- dim(catch)[1]
    catch_sets <- dim(catch)[2]
    }else{
      len_catch <- length(catch)
      catch_sets <- 1
    }

  res_months <- data.frame(month = 1:12)
  # number of recruits per month
  N_months_loop <- matrix(NA, nrow = 12, ncol = catch_sets)
  ti <- matrix(NA, nrow = len_catch, ncol = catch_sets)
  tS_frac <- matrix(NA, nrow = len_catch, ncol = catch_sets)
  correspond_month <- matrix(NA, nrow = len_catch, ncol = catch_sets)
  for(i in 1:catch_sets){
    if(catch_sets > 1){
      catchi <- catch[,i]
      tsampli <- tsample[i]
    }else{
      catchi <- catch
      tsampli <- tsample
    }

    # # special vBGF (D = 1)
    # ti[,i] <- log(1 - classes.num/Linf)/-K + t0 - tsampli
    # # generalised vBGF
    # if(!is.na(D)) ti[,i] <- log(1 - (classes.num/Linf)^D)/(-K*D) + t0 - tsampli
    # # seasonalized vBGF
    # if(!is.na(C) & !is.na(ts)) ti[,i] <- (log(1 - (classes.num/Linf)^D) +
    #                                     C * ((K*D)/2*pi)*sin(2*pi)*((t0-tsampli)-ts))/(-K*D) + t0 - tsampli

    ti[,i] <- (log(1 - (classes.num/Linf)^D) +
                 C * ((K*D)/2*pi)*sin(2*pi)*((t0-tsampli)-ts))/(-K*D) + t0 - tsampli

    # t at S = 0 as fraction of year
    tS_frac[,i] <- ti[,i] - floor(ti[,i])
    # corresponding months
    correspond_month[,i] <- floor(tS_frac[,i] * 12 + 1)

    N_months_pre <- aggregate(list(numbers=catchi), by = list(month = correspond_month[,i]),
                              sum, na.rm =TRUE)
    N_months_all <- merge(res_months, N_months_pre, by.x = "month", all.x=TRUE)
    N_months_all$numbers[is.na(N_months_all$numbers)] <- 0
    N_months_loop[,i] <- N_months_all$numbers
  }

  N_months_final <- apply(N_months_loop, 1, FUN = mean, na.rm = TRUE)

  # percentage of recruits per month
  N_months_per <- (N_months_final/sum(N_months_final, na.rm = TRUE)) * 100

  if("t0" %in% names(res)){
    months_abb <- month.abb[res_months$month]
  }

  if(catch_sets == 1){
    ti <- ti[,1]
    tS_frac <- tS_frac[,1]
    correspond_month <- correspond_month[,1]
    N_months_loop <- N_months_loop[,1]
  }

  res2 <- c(res,list(ti = ti,
                    tS_frac = tS_frac,
                    cor_months = correspond_month,
                    months = res_months$month))
  if("t0" %in% names(res)) res2$months_abb <- months_abb
  ret <- c(res2,list(all_recruits = N_months_loop,
                    mean_recruits = N_months_final,
                    per_recruits = N_months_per))

  class(ret) <- "recruitment"
  if(plot) plot(ret)
  return(ret)
}
