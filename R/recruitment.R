#' @title Recruitment patterns
#'
#' @description  This function estimates recrutiment patterns from length-frequency data-
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \code{midLengths}: midpoints of the length classes (length-frequency
#'   data),
#'   \item \code{Linf}: infinite length for investigated species in cm [cm],
#'   \item \code{K}: growth coefficent for investigated species per year [1/year],
#'   \item \code{t0}: theoretical time zero, at which individuals of this species hatch,
#'   \item \code{D}: optional; for generalised vBGF
#'   \item \code{C}: optional; intensity of the (sinusoid) growth oscillations of
#'   the model,
#'   \item \code{ts}: optional; onset of the first oscillation with regard to t = 0,
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
#' data(trout)
#'
#' # add growth parameters
#' trout$Linf <- 16
#' trout$K <- 0.4
#'
#' # retrieve sampling times from catch matrix
#' s_dates <- strsplit(colnames(trout$catch),"X")
#' s_dates <- unlist(lapply(s_dates, function(x) return(x[2])))
#' s_dates <- as.POSIXlt(s_dates, format="%d.%m.%Y")
#'
#' recruitment(trout, tsample = s_dates$yday/365)
#'
#' @details
#' This function allows to extract information about the recruitment patterns of a
#' stock by backward projection onto the length axis of a set of length frequency data
#' using the (special, generalised or seasonalised)
#' von Bertallanfy growth curve (vBGF, Pauly 1982). The method assumes that (i) all fish in a
#' data set grow as described by a single set of growth parameters and (ii) one month out
#' of twelve always has zero recruitment. The second assumption is probably not met as
#' in temperate waters more than one month has zero recruitment and in the tropics
#' there might always be recruitment.
#' If t0 is not provided, only a realtive recruitment pattern can be estimated. It is
#' therefore not possible to assign months to the results. However, an estimate of t0
#' can be obtained by the time
#' lag between peak spawning and recruitment. Several length-frequency data sets
#' can be used to estimate the recrutiment pattern by providing catch as a matrix and
#' setting catch_column to NA (default). Then the fraction per time is calculated for
#' each size class in each sample an then pooled together. For the generalised vBGF D is
#' required, for the seasonalised vBGF C, ts and D.
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{ti}: actual age,
#'   \item \strong{tS_frac}: age at which the lenght was zero expressed as fraction
#'   of the year,
#'   \item \strong{correspond_month}: corresponding month,
#'   \item \strong{N_months_all}: dataframe with months and the number of recruits
#'   per month,
#'   \item \strong{N_months_per}: precentage of number of recruits per month.
#' }
#'
#' @references
#' Brey, T., Soriano, M., Pauly, D., 1988. Electronic length frequency analysis. A
#' revised and expanded user's guide to ELEFAN 0, 1 and 2. (Second edition).
#' Berichte aus dem Institut für Meereskunde Kiel, No 177, 31p.
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
  D <- ifelse("D" %in% names(res),res$D, NA)
  C <- ifelse("C" %in% names(res),res$C, NA)
  ts <- ifelse("ts" %in% names(res),res$ts, NA)
  classes <- as.character(res$midLengths)
  # create column without plus group (sign) if present
  classes.num <- do.call(rbind,strsplit(classes, split="\\+"))
  classes.num <- as.numeric(classes.num[,1])

  if(is.na(catch_column)) catch <- res$catch
  if(!is.na(catch_column)) catch <- res$catch[,catch_column]

  # special vBGF (D = 1)
  ti <- log(1 - classes.num/Linf)/-K + t0 - tsample
  # generalised vBGF
  if(!is.na(D)) ti <- log(1 - (classes.num/Linf)^D)/(-K*D) + t0 - tsample
  # seasonalized vBGF
  if(!is.na(C) & !is.na(ts)) ti <- (log(1 - (classes.num/Linf)^D) +
                                      C * ((K*D)/2*pi)*sin(2*pi)*((t0-tsample)-ts))/(-K*D) + t0 - tsample

  # t at S = 0 as fraction of year
  tS_frac <- ti - floor(ti)
  # corresponding months
  correspond_month <- floor(tS_frac * 12 + 1)


  # results
  res_months <- data.frame(month = 1:12)
  # number of recruits per month
  N_months <- aggregate(list(numbers=catch), by = list(month = correspond_month),
                        sum, na.rm =TRUE)

  N_months_all <- merge(res_months, N_months, by.x = "month", all.x=TRUE)
  N_months_all$numbers[is.na(N_months_all$numbers)] <- 0

  # percentage of recruits per month
  N_months_per <- (N_months_all$numbers/sum(N_months_all$numbers, na.rm = TRUE)) * 100


  ret <- c(res,list(ti = ti,
                    tS_frac = tS_frac,
                    correspond_month = correspond_month,
                    N_months_all = N_months_all,
                    N_months_per = N_months_per))
  class(ret) <- "recruitment"
  if(plot) plot(ret)
  return(ret)
}
