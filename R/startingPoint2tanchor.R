#' @title Convert FiSAT's starting point to t_anchor value
#'
#' @description Starting points returned or chosen within FiSAT are not supported
#'    in TropFishR. Instead \code{t_anchor} takes on the job of anchoring VBGF growth curves on
#'    a temporal axis. This function allows to convert FiSAT's starting points to \code{t_anchor} values
#'
#' @param param lfq object with dates, midLengths, catch, and growth parameters 'Linf'
#'    and 'K' from VBGF
#' @param startingLength starting length as returned by FiSAT, indicating the length within
#'    the starting sample cut by a growth curve
#' @param startingSample starting sample as returned by FiSAT, indicating the sample which is
#'    cut by a growth curve
#'
#' @keywords function lfq startingPoints t_anchor
#'
#' @examples
#' data(synLFQ1)
#' synLFQ1$Linf <- 54
#' synLFQ1$K <- 0.6
#' lfqNEW <- startingPoint2tanchor(synLFQ1, startingLength = 20, startingSample = 4)
#'
#' @return lfq object with t_anchor value
#'
#' @export

startingPoint2tanchor <- function(param, startingLength, startingSample){

  res <- param
  if(("Linf" %in% names(res)) == FALSE) stop(noquote("Please add Linf to param!"))
  if(("K" %in% names(res)) == FALSE) stop(noquote("Please add K to param!"))
  if(("dates" %in% names(res)) == FALSE) stop(noquote("Please add a dates to param!"))

  if(length(res$dates) < startingSample) stop(noquote(paste("Only",
                                                            length(res$dates),
                                                            "sample dates in param$dates, but startingSample is",
                                                            startingSample)))

  tx <- res$dates[startingSample]
  tx <- date2yeardec(tx)
  agex <- VBGF(L = startingLength, param = list(Linf = res$Linf, K = res$K, t0 = 0))
  res$t_anchor <- (tx-agex) %% floor(tx-agex)

  return(res)
}
