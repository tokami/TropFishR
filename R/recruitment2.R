

#' Calculate recruitment pattern of lfq object
#'
#' \code{recruitment2} calculates the recruitment pattern (i.e. monthly)
#' of an lfq object \code{\link{VBGF}} parameters and the length at
#' recruitment (\code{Lrecr}, e.g. to the population or the fishery).
#' Time at recruitment \code{trecr} is calculated via a hindcast
#' using the \code{\link{calc_trecr}} function.
#'
#' This function does not represent a robust statistical method as it
#' assumes the same growth parameters for all individual counts
#' in the lfq count bins.
#' Nevertheless, the results should provide general information about
#' the recruitment pattern, e.g. relative recruitment strength by month.
#'
#' @param lfq a length-frequency object (i.e. class \code{lfq})
#' @param Lrecr numeric. Length at recruitment
#' (default: \code{Lrecr = 0}).
#' @param plot logical. Should monthly recruitment pattern by plotted
#' (default: \code{plot = TRUE}).
#' @param hide.progressbar logical. Should progress bar be displayed
#' (default: \code{hide.progressbar = FALSE})
#'
#' @return object of "histogram" class containing
#' monthly recruitment pattern
#' @export
#'
#' @examples
#' lfq <- synLFQ4
#' lfq$par <- list(
#'   Linf = 80, K = 0.5, t_anchor = 0.25, C = 0.75, ts = 0.5
#' )
#' tmp <- recruitment2(lfq = lfq)
#' tmp
#'
recruitment2 <- function(
  lfq,
  Lrecr = 0,
  plot = TRUE,
  hide.progressbar = FALSE
){
  trecr <- NaN * lfq$catch
  Lincl <- which(lfq$midLengths < lfq$par$Linf & lfq$midLengths > Lrecr)
  trecr <- trecr[Lincl,]
  counter_mat <- matrix(seq(length(trecr)), nrow = nrow(trecr), ncol = ncol(trecr), byrow = TRUE)
  if(!hide.progressbar){
    pb <- txtProgressBar(min=1, max=length(trecr), style=3)
  }
  for(i in seq(nrow(trecr))){
    for(j in seq(ncol(trecr))){
      trecr[i,j] <- calc_trecr(
        par = lfq$par,
        Lrecr = Lrecr,
        Lt = lfq$midLengths[i],
        t = lfq$dates[j]
      )
      setTxtProgressBar(pb, counter_mat[i,j])
    }
  }
  if(!hide.progressbar){close(pb)}
  trecr <- rep(x = c(trecr), times=c(lfq$catch[Lincl,]))
  trecr <- trecr[!is.na(trecr)]
  trecr <- yeardec2date(trecr)
  trecr.mo <- as.numeric(format(trecr, "%m"))
  h <- hist(trecr.mo, breaks = seq(0.5, 12.5, by=1), plot = FALSE)
  if(plot){
    barplot(h$counts/sum(h$counts)*100, names.arg = month.abb, las = 2, ylab = "Recruitment [%]")
  }
  return(h)
}
