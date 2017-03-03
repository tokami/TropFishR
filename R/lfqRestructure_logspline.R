#' @title Restructuring of length frequency data using log-transformed frequencies
#'
#' @description This is an alternate (experimantal) LFQ restructuring method
#' that works on log-transformed frequencies, followed by a calculation
#' of deviations againsta a fitted spline function.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes
#'   \item \strong{dates} dates of sampling times (class Date)
#'   \item \strong{catch} matrix with catches/counts per length class (row) and
#'   sampling date (column)
#' }
#' @param n.mult multiplier to use in expanding the number of bins (Default=5)
#' @param n.peak expected number of cohorts (Default=5)
#' @param endrule running median rule (Default="median")
#' @param logscore Logical. Return log-transformed scores (Default=FALSE)
#'
#' @examples
#' data("synLFQ4")
#'
#' synLFQ4 <- lfqRestructure_logspline(synLFQ4, n.peak=5)
#' plot(synLFQ4)
#'
#' synLFQ4 <- lfqRestructure_logspline(synLFQ4, n.peak=7, logscore=TRUE)
#' plot(synLFQ4)
#'
#' # compare to lfqRestructure
#' op <- par(mfcol=c(2,1), mar=c(4,4,1,1))
#' synLFQ4 <- lfqRestructure(synLFQ4, MA=13)
#' plot(synLFQ4)
#' synLFQ4 <- lfqRestructure_logspline(synLFQ4, n.peak=7, logscore=FALSE)
#' plot(synLFQ4)
#' par(op)
#'
#' @export

lfqRestructure_logspline <- function(
  param,
  n.mult = 5,
  n.peak = 5,
  endrule = "median",
  logscore = FALSE
){

  nearest <- function (value, lookup_vector){
  	diff_sq=(value-lookup_vector)^2
  	order(diff_sq, decreasing = FALSE)[1]
  }

  lfq <- param

  x <- lfq$midLengths
  logx <- log(x)
  rcounts <- NaN*lfq$catch
  for(i in seq(ncol(lfq$catch))){
    y <- lfq$catch[,i]
    k <- length(x) * n.mult / n.peak
    spl <- spline(logx, y, n = length(x)*n.mult)
    runm <- suppressWarnings(runmed(spl$y, k = k, endrule = endrule))
    runm <- runm[unlist(lapply(as.list(logx), FUN = function(x){nearest(x, spl$x)}))]
    score <- y - runm
    if(logscore){
      score <- (log(score - min(score) + 1) - log(0 - min(score) + 1))
    }
    rcounts[,i] <- score
  }
  lfq$rcounts <- rcounts

  # create peak matrix
  prep_mat <- lfq$rcounts
  prep_mat <- ifelse(prep_mat > 0,1,0)
  peaks_mat <- NA*prep_mat
  for(i in seq(ncol(peaks_mat))){
    vec_peaki <- prep_mat[,i]
    runs <- rle(vec_peaki)
    rle_val <- runs$values
    rle_val[which(rle_val == 1)] <- 1:length(rle_val[which(rle_val == 1)])
    peaks_mat[,i] <- rep(rle_val, runs$lengths)
  }
  maxn.peaks <- max(peaks_mat, na.rm=TRUE)
  peaks_mat <- peaks_mat + (prep_mat * maxn.peaks * col(peaks_mat))
  lfq$peaks_mat <- peaks_mat

  # ASP calc
  sampASP <- NaN*seq(ncol(rcounts))
  for(i in seq(ncol(rcounts))){
    # lfq.i <- lfq[i,]
    tmp <- rle(sign(rcounts[,i]))
    start.idx <- c(1, cumsum(tmp$lengths[-length(tmp$lengths)])+1)
    end.idx <- cumsum(tmp$lengths)
    posrun <- which(tmp$values == 1)
    peakval <- NaN*posrun
    for(p in seq(length(posrun))){
      peakval[p] <- max(rcounts[start.idx[posrun[p]]:end.idx[posrun[p]], i ])
    }
    sampASP[i] <- sum(peakval)
  }
  ASP <- sum(sampASP)
  lfq$ASP <- ASP
  lfq$MA <- NA

  class(lfq) <- "lfq"
  return(lfq)
}
